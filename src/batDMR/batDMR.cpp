#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <string.h>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iterator>
#include "stddef.h"
#include <gsl/gsl_cdf.h>
#include <tr1/cmath>
#include <gsl/gsl_sf_gamma.h>
#include "regression.hpp"
#include "combine_pvals.hpp"
#include "merge.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include "zlib.h"
#if __GNUC__>2
//#include <ext/hash_set>
#include <ext/hash_map>
using namespace __gnu_cxx;
#else
//#include <hash_set>
#include <hash_map>
using namespace stdext;
#endif


#define BUFSIZE 1000000

using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::ios_base;
using std::cout;
using std::ifstream; 
using std::ofstream; 
using std::floor;
using std::map;
map <string,int> String_Hash;
bool geneid=false;
int Genome_Count = 0;
std::string filtercontext = "ALL";

struct str_hash{
        size_t operator()(const string& str) const
        {
                return __stl_hash_string(str.c_str());
        }
};

struct compare_str{
        bool operator()(const std::string& p1, const std::string& p2) const{
                return p1 == p2;
        }
};

typedef hash_map<std::string, int, str_hash, compare_str> chrposMap;

/*
 *this paragram part from cufdiff: Benjamani-Hochberg procedure and fisher test studied bedtools. in addation, the beta-binomail distribution from RADmeth software 
*/

void chromLengthExact(string & refSeqFile)
{
        FILE* iFPtr;
        FILE* oFPtr;
        
        string chrLenFile;
        chrLenFile=refSeqFile;
        chrLenFile+=".len";
        
        fprintf(stderr, "refSeqFile = %s\n", refSeqFile.c_str());
        fprintf(stderr, "chrLenFile = %s\n", chrLenFile.c_str());
        
        char* readBuffer = (char*)malloc(sizeof(char) * BUFSIZE);
        int len = 0;
        int lines = 0;
        
        char* token;
        char seps[] = " \t\n\r";
        char chrName[100];
        
        if(!(iFPtr = fopen(refSeqFile.c_str(), "r"))) {
                printf("File %s open error.\n", refSeqFile.c_str());
                exit(1);
        }
        
        if(!(oFPtr = fopen(chrLenFile.c_str(), "r"))) {
                oFPtr = fopen(chrLenFile.c_str(), "w");
                while(fgets(readBuffer, BUFSIZE, iFPtr)) {
                if(strlen(readBuffer) >= BUFSIZE - 1) {
                           fprintf(stderr, "Too many characters in one row! Try to split the long row into several short rows (fewer than %d characters per row).\n", BUFSIZE);
                           exit(1);
                }

                if(readBuffer[0] == '>') {
                        if(lines > 0) {
                                fprintf(oFPtr, "%s\t%d\n", chrName, len);
                        }
                        // Save name
                        token = strtok(readBuffer + 1, seps);
                        strcpy(chrName, token);
                        len = 0;
                }
                else {
                        // Substract \n
                        len += strlen(readBuffer) - 1;
                }

                lines++;
                }

             if(lines > 0) {
                fprintf(oFPtr, "%s\t%d\n", chrName, len);
             }
            if(!(oFPtr = fopen(chrLenFile.c_str(), "r"))){
                printf("File %s open error.\n", chrLenFile.c_str());
                exit(1);
            }
        }
        
        fclose(iFPtr);
        fclose(oFPtr);
        
        return;
};

void Print_dm_result(const vector<PvalLocus> &pval_loci, ostream &output_encoding, double cutoff, double Pcutoff, double methdiff, string dmr_outfile, 
    bool singleAuto, int mindmc, int mindmcdis, int maxdmrlen) {
  string record, chrom, context, sign;string name;
  size_t position, coverage_factor, meth_factor, coverage_rest, meth_rest;
  double pval;double adjust_pvalue;
  ofstream OutFiledmr;
  if(singleAuto) OutFiledmr.open(dmr_outfile.c_str());

  vector<PvalLocus>::const_iterator cur_locus_iter = pval_loci.begin();
  int Ndmc = 0; int prevdmc = 0; int dmrtotalLen = 0; string prevchrom = "null";
  int dmrC1 = 0; int dmrCover1 = 0; int dmrC2 = 0; int dmrCover2 = 0;
  int dmrstart = 0; int dmrend = 0;
  int hyperNdmc=0,hypoNdmc=0;
  for(unsigned i=0;i<pval_loci.size();) 
  {
    pval= cur_locus_iter->raw_pval;
    adjust_pvalue=cur_locus_iter->adjust_pval; 
    coverage_factor=cur_locus_iter->coverage_factor;
    meth_factor=cur_locus_iter->meth_factor;
    coverage_rest=cur_locus_iter->coverage_rest;
    meth_rest=cur_locus_iter->meth_rest;
    position = cur_locus_iter->position;
    chrom = cur_locus_iter->chrom;
    sign = cur_locus_iter-> sign ; context = cur_locus_iter->context ; //std::string state="Faild";
    std::string name=cur_locus_iter->name.c_str();
    double meth_diff=fabs(double(meth_factor)/coverage_factor - double(meth_rest)/coverage_rest) ;
    double signed_methdiff=double(meth_factor)/coverage_factor - double(meth_rest)/coverage_rest;
    //if(!(adjust_pvalue < cutoff && (pval < Pcutoff || (pval < Pcutoff+0.05 && coverage_factor+coverage_rest<=50) ) && meth_diff >= methdiff ) ) {
    if(!(adjust_pvalue < cutoff && pval < Pcutoff && meth_diff >= methdiff ) ) {
        cur_locus_iter++;i++;
        continue;
    };
    
    output_encoding << chrom << "\t" << position << "\t" << sign << "\t"
                    << context << "\t" << pval << "\t";
    if (0 <= pval && pval <= 1) {
      output_encoding << adjust_pvalue << "\t";
      cur_locus_iter++;i++;
    } else {
      output_encoding << pval << "\t";
      cur_locus_iter++;i++;
    }
    
    //if(geneid) name=cur_locus_iter->name.c_str();
    if(geneid) output_encoding << meth_factor  << "\t" << coverage_factor << "\t"
                    << meth_rest  << "\t" << coverage_rest << "\t" << signed_methdiff << "\t" <<  name.c_str() << "\n" ;
    else {
        output_encoding << meth_factor  << "\t" << coverage_factor << "\t"<< meth_rest  << "\t" << coverage_rest << "\t" << signed_methdiff
    			<< "\n" ;
        if(singleAuto){
            if(prevchrom=="null") prevchrom = chrom;
//fprintf(stderr, "\nHHHH %d %d %d %d %d %d %d", Ndmc, mindmc, dmrstart, dmrend, prevdmc, position, position - prevdmc );
            if(dmrstart == 0){ 
                dmrstart = position;
                prevdmc = position;
                Ndmc=1 ;
                if(signed_methdiff > 0) {
                     hyperNdmc=1;
                     hypoNdmc=0;
                }
                else{
                     hypoNdmc=1;
                     hyperNdmc=0;
                }
            }
            else if(prevchrom != chrom) {
//fprintf(stderr, "\nzzzTTT %d %d %d %s %s", dmrtotalLen, dmrtotalLen + position - prevdmc, maxdmrlen, prevchrom.c_str(), chrom.c_str());
                if(Ndmc >= mindmc){
                    OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                    << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                    << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                }
                Ndmc = 0; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                dmrC1 = 0; dmrCover1 = 0; dmrC2 = 0; dmrCover2 = 0;
                dmrstart = position;
                prevchrom = chrom;
                hyperNdmc = 0; hypoNdmc = 0;
            }
            else if(position - prevdmc <= mindmcdis) {
//fprintf(stderr, "\nTTT %d %d %d", dmrtotalLen, dmrtotalLen + position - prevdmc, maxdmrlen);
               if(dmrtotalLen + position - prevdmc <= maxdmrlen){
                   Ndmc ++ ;
                   if(signed_methdiff > 0) hyperNdmc++;
                   else hypoNdmc++;
                   dmrend = position;
                   dmrtotalLen = dmrtotalLen + position - prevdmc;
                   prevdmc = position;
                   dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                   dmrC2 += meth_rest; dmrCover2 += coverage_rest;
               }else{
                   if(Ndmc >= mindmc){
                       OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                        << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                        << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                   }
                   Ndmc = 1; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                   if(signed_methdiff > 0) {
                        hyperNdmc=1;
                        hypoNdmc=0;
                   }
                   else{
                        hypoNdmc=1;
                        hyperNdmc=0;
                   }
                   dmrstart = position;
                   dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                   dmrC2 += meth_rest; dmrCover2 += coverage_rest;
               }
            }else{ // jian ge tai da, should be here
                if(Ndmc >= mindmc){
                     OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
                     << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
                     << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
                }
                Ndmc = 1; dmrend = 0; prevdmc = position; dmrtotalLen = 0;
                if(signed_methdiff > 0) {
                     hyperNdmc=1;
                     hypoNdmc=0;
                }
                else{
                     hypoNdmc=1;
                     hyperNdmc=0;
                }
                dmrstart = position;
                dmrC1 += meth_factor; dmrCover1 += coverage_factor;
                dmrC2 += meth_rest; dmrCover2 += coverage_rest;
            }
        }// end single auto
    } //end print
  } //end for
//printf("\n%d\n", Ndmc);
  if(Ndmc >= mindmc){
     OutFiledmr << prevchrom << "\t" << dmrstart << "\t" << dmrend << "\t"
     << (double)dmrC1/dmrCover1 << "\t" << (double)dmrC2/dmrCover2 << "\t" << Ndmc
     << "\t" << hyperNdmc << "," << hypoNdmc << "\n";
  }
  return;
};

/*
bool p_value_lt(const SampleDifference* lhs, const SampleDifference* rhs)
{
	return lhs->p_value < rhs->p_value;
}

// Benjamani-Hochberg procedure
int fdr_significance(double fdr, 
                     vector<SampleDifference*>& tests)
{
	sort(tests.begin(), tests.end(), p_value_lt);
	vector<SampleDifference*> passing;
    
	for (int k = 0; k < (int)tests.size(); ++k)
	{
		if (tests[k]->test_status == OK)
		{
			passing.push_back(tests[k]);
		}
		else
		{
			tests[k]->significant = false;
		}
	}
	int significant = 0;
	float pmin=1;
	int n = (int) passing.size();
    //use the same procedure as p.adjust(...,"BH") in R
	for (int k = n-1; k >= 0; k--)
	{
		double corrected_p = (double) passing[k]->p_value * ((double) n/(double) (k+1));
        //make sure that no entry with lower p-value will get higher q-value than any entry with higher p-value
		if (corrected_p < pmin) 
		{
			pmin = corrected_p;
		}
		else
		{
			corrected_p = pmin;
		}
        // make sure that the q-value is always <= 1 
		passing[k]->corrected_p = (corrected_p < 1 ? corrected_p : 1); 
		passing[k]->significant = (corrected_p <= fdr);
        significant += passing[k]->significant;
	}
    
	return passing.size();
}
*/
static inline double log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
};

// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
          gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
          (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
};


static double fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
};

struct DNAmethException {
  DNAmethException(std::string m) : message(m) {}
  std::string what() const {return message;}
  std::string message;
};

class Cpg {
public:
  Cpg() {};
  Cpg(std::string encoding);
  void set_chrom(std::string chrom) { chrom_ = chrom; }
  std::string chrom() const { return chrom_; }
  void set_locus(size_t locus) { locus_ = locus; }
  size_t locus() const { return locus_; }
  void set_strand(std::string strand) { strand_ = strand; }
  std::string strand() const { return strand_; }
  void set_context(std::string context) { context_ = context; }
  std::string context() const { return context_; }
  void set_meth(size_t meth) { meth_ = meth; }
  size_t meth() const { return meth_; }
  void set_total(size_t total) { total_ = total; }
  size_t total() const { return total_; }
  void set_name(std::string name) { name_ = name; }
  std::string name() const { return name_; }
private:
  std::string chrom_;
  size_t locus_;
  std::string strand_;
  std::string context_;
  size_t total_;
  size_t meth_;
  std::string name_;
};

Cpg::Cpg(string encoding) {
  //char sign;
  //string name;
  //double level;
    istringstream encoding_stream(encoding);	
	//if(!(encoding_stream >> chrom_ >> locus_ >> strand_ >> context_ >> meth_ >> total_ >> name))
	if(encoding[0]=='#') return;
	  if (!(encoding_stream >> chrom_ >> locus_ >> strand_ >> context_ >> meth_ >> total_))
		throw (std::logic_error("Couldn't parse a line \"" + encoding + "\".\n"));
	if(geneid) encoding_stream >> name_;
};

struct Methy_Hash
{
	//int *positiveSample;
	//int *Neg_positiveSample;
	chrposMap positiveMap;
	chrposMap Neg_positiveMap;
	int Index;
};

std::string int_to_String(int n)
{
    int max = 100;
    int m = n;
    char s[max];
    char ss[max];
    int i=0,j=0;
    if (n < 0)
    {
        m = 0 - m;
        j = 1;
        ss[0] = '-';
    }    
    while (m>0)
    {
        s[i++] = m % 10 + '0';
        m /= 10;
    }
    s[i] = '\0';
    i = i - 1;
    while (i >= 0)
    {
        ss[j++] = s[i--];
    }    
    ss[j] = '\0';    
    return ss;
}

void merge_methylomes_gz(vector<string> names, vector<string> methylomes, ostream &count_table, ostream &count_inform_table , unsigned Sample1Size,Methy_Hash* Methy_List,std::map <string,int> &String_Hash) {
  
  vector<string>::const_iterator it = names.begin();
  
  count_inform_table << "base\tcase\n" << *it;
  count_table << *it++;
  unsigned nSize=0;
  while (it != names.end())
  {
  if(nSize < Sample1Size)
  { 
    count_inform_table << "\t1\t0\n" << *it;
    nSize++;
  }else
    count_inform_table << "\t1\t1\n" << *it;
  count_table << "\t" << *it++;
  }
  count_inform_table << "\t1\t1" ;
  
  //
  int MAXDES = 600;
  //output matrix
  vector<gzFile> methylomes_fstream;
  vector<string>::iterator meth_it = methylomes.begin();
  while(meth_it != methylomes.end())
  {
      methylomes_fstream.push_back(gzopen((*meth_it).c_str(), "rb"));
      ++meth_it;
  }

  //vector<gzFile>::iterator meth_if = methylomes_fstream.begin();
  char encoding[MAXDES];
  for(int i=0;i<methylomes_fstream.size();i++){
    strcpy(encoding, "0");
    gzgets(methylomes_fstream[i],encoding,MAXDES);
    if(encoding[0]=='0' || strlen(encoding)<5){
      fprintf(stderr, "\nnot find input file!\n");
      exit(0);
    }
    //if(meth_if == methylomes_fstream.end() ) break;
    //++meth_if;
  }

  std::string methTmp;
  while(true)
  {
    vector<gzFile>::iterator meth_it = methylomes_fstream.begin();
  char encoding[MAXDES];
  while(meth_it != methylomes_fstream.end() ) {
    if (!gzgets(*meth_it,encoding,MAXDES)) 
    {
      meth_it = methylomes_fstream.erase(meth_it);
      continue;
    }
    if(encoding[0]=='#') continue;
    Cpg cpg((string)encoding);

    if(String_Hash.find(cpg.chrom())==String_Hash.end() || (cpg.context() != filtercontext && filtercontext != "ALL") ) 
    {
      //meth_it++;
      continue;
    }

    int H=String_Hash[cpg.chrom()];
    if(cpg.total()>0) 
    {
      if(!geneid){
        methTmp = cpg.chrom() + "_";
        methTmp = methTmp + int_to_String(cpg.locus());
      }else methTmp = cpg.name();
      if(cpg.strand()=="+"){
        Methy_List[H].positiveMap[methTmp]++; 
      }else if(cpg.strand()=="-"){
        Methy_List[H].Neg_positiveMap[methTmp]++;
      }
    }
    //meth_it++;
  }

    if( methylomes_fstream.size()==0 )
    {
      methylomes_fstream.clear();
      break;
    }
  }

  meth_it = methylomes.begin();
  while(meth_it != methylomes.end())
  {
      methylomes_fstream.push_back(gzopen((*meth_it).c_str(), "rb"));
      meth_it++;
  }

  while (true) {
    
    vector<gzFile>::iterator meth_it = methylomes_fstream.begin();
    char encoding[MAXDES];
    if (!gzgets(*meth_it,encoding,MAXDES))
      break;
      
    Cpg cpg1((string)encoding);
    if(String_Hash.find(cpg1.chrom())==String_Hash.end()  || (cpg1.context() != filtercontext && filtercontext != "ALL")) continue;
    count_table << "\n";
    int H=String_Hash[cpg1.chrom()];
    unsigned int loci=cpg1.locus();
    string strand=cpg1.strand();
    if(!geneid){
      methTmp = cpg1.chrom() + "_";
      methTmp = methTmp + int_to_String(loci);
    }else methTmp = cpg1.name();
    while(  (strand=="+" && Methy_List[H].positiveMap[methTmp] != (int)methylomes.size())  || 
        (strand=="-" && Methy_List[H].Neg_positiveMap[methTmp] != (int)methylomes.size() )  )
    {
      if (!gzgets(*meth_it,encoding,MAXDES))
        break;
      Cpg cpg((string)encoding);
      cpg1=cpg;
      if(String_Hash.find(cpg1.chrom())==String_Hash.end()  || (cpg1.context() != filtercontext && filtercontext != "ALL")) continue;
  if(!geneid){
    methTmp = cpg1.chrom() + "_";
    methTmp = methTmp + int_to_String(cpg1.locus());
  }else methTmp = cpg1.name();
      H=String_Hash[cpg1.chrom()];
      loci=cpg1.locus();
      strand=cpg1.strand();
    }
    
    bool validc = false;    
    if(  ( (strand=="+" && Methy_List[H].positiveMap[methTmp] == (int)methylomes.size())  || 
        (strand=="-" && Methy_List[H].Neg_positiveMap[methTmp] == (int)methylomes.size() ) ) && cpg1.total()>0) {
      count_table << cpg1.chrom() << ":" << cpg1.locus() << ":" << cpg1.strand() << ":" << cpg1.context() 
                << "\t" << cpg1.meth() << "\t" << cpg1.total();
        validc=true;
    }

    meth_it++;
    while(validc && meth_it != methylomes_fstream.end()) {

  while( gzgets(*meth_it,encoding,MAXDES) )
  {
          Cpg cpg2((string)encoding);
          if( (cpg2.context() != filtercontext && filtercontext != "ALL")) continue;
          if( (geneid && cpg2.name()==cpg1.name()) || (!geneid && H==String_Hash[cpg2.chrom()] &&  cpg1.locus()==cpg2.locus() && cpg1.strand()==cpg2.strand() && cpg2.total()>0) )
          {
          count_table << "\t" << cpg2.meth() << "\t" << cpg2.total();
          if(geneid) count_table << "\t" << cpg2.name();
          break;
          }
    }
      meth_it++;
    }
  }

  for(size_t ind = 0; ind < methylomes_fstream.size(); ++ind)
    delete methylomes_fstream[ind];
}
//
void merge_methylomes(vector<string> names, vector<string> methylomes, ostream &count_table, ostream &count_inform_table , unsigned Sample1Size,Methy_Hash* Methy_List,std::map <string,int> &String_Hash) {
  
  vector<string>::const_iterator it = names.begin();
  
  count_inform_table << "base\tcase\n" << *it;
  count_table << *it++;
  unsigned nSize=0;
  while (it != names.end())
  {
	if(nSize < Sample1Size)
	{	
		count_inform_table << "\t1\t0\n" << *it;
		nSize++;
	}else
		count_inform_table << "\t1\t1\n" << *it;
	count_table << "\t" << *it++;
  }
  count_inform_table << "\t1\t1" ;
  
  //output matrix
  vector<istream*> methylomes_fstream;
  vector<string>::iterator meth_it = methylomes.begin();
  while(meth_it != methylomes.end())
  {
  	  methylomes_fstream.push_back(new ifstream( (*meth_it).c_str() ));
  	  ++meth_it;
  }

  vector<istream*>::iterator meth_if = methylomes_fstream.begin();
  while(meth_if != methylomes_fstream.end()){
    string encoding;
    getline(*(*meth_if), encoding);
    if (encoding.empty())
    {
        fprintf(stderr, "no such file, please check the file name!\n");
        exit(0);
    }
    meth_if++;
  }

  std::string methTmp;
  while(true)
  {
  	vector<istream*>::iterator meth_it = methylomes_fstream.begin();
	string encoding;
	while(meth_it != methylomes_fstream.end() ) {
		getline(*(*meth_it), encoding);
		if (encoding.empty()) 
		{
			meth_it = methylomes_fstream.erase(meth_it);
			continue;
		}
                if(encoding[0]=='#') continue;
		Cpg cpg(encoding);
		if(String_Hash.find(cpg.chrom())==String_Hash.end()  || (cpg.context() != filtercontext && filtercontext != "ALL")) 
		{
			//meth_it++;
			continue;
		}
		int H=String_Hash[cpg.chrom()];
		if(cpg.total()>0) 
		{
			if(!geneid){
				methTmp = cpg.chrom() + "_";
				methTmp = methTmp + int_to_String(cpg.locus());
			}else methTmp = cpg.name();
			if(cpg.strand()=="+"){
				Methy_List[H].positiveMap[methTmp]++;	
			}else if(cpg.strand()=="-"){
                                Methy_List[H].Neg_positiveMap[methTmp]++;
			}
		}
		//meth_it++;
	}

	if( methylomes_fstream.size()==0 )
	{
		methylomes_fstream.clear();
		break;
	}
  }
  meth_it = methylomes.begin();
  while(meth_it != methylomes.end())
  {
  	  methylomes_fstream.push_back(new ifstream( (*meth_it).c_str() ));
  	  meth_it++;
  }
  
  while (true) {
    
    vector<istream*>::iterator meth_it = methylomes_fstream.begin();
    string encoding;
    getline(*(*meth_it), encoding);
    if (encoding.empty())
      break;
      
    Cpg cpg1(encoding);
    if(String_Hash.find(cpg1.chrom())==String_Hash.end()  || (cpg1.context() != filtercontext && filtercontext != "ALL")) continue;
    count_table << "\n";
    int H=String_Hash[cpg1.chrom()];
    unsigned int loci=cpg1.locus();
    string strand=cpg1.strand();
    if(!geneid){
	    methTmp = cpg1.chrom() + "_";
	    methTmp = methTmp + int_to_String(loci);
    }else methTmp = cpg1.name();
    while(  (strand=="+" && Methy_List[H].positiveMap[methTmp] != (int)methylomes.size())  || 
    		(strand=="-" && Methy_List[H].Neg_positiveMap[methTmp] != (int)methylomes.size() )  )
    {
    	getline(*(*meth_it), encoding);
    	if (encoding.empty())
      	break;
    	Cpg cpg(encoding);
    	cpg1=cpg;
    	if(String_Hash.find(cpg1.chrom())==String_Hash.end()  || (cpg1.context() != filtercontext && filtercontext != "ALL")) continue;
	if(!geneid){
		methTmp = cpg1.chrom() + "_";
		methTmp = methTmp + int_to_String(cpg1.locus());
	}else methTmp = cpg1.name();
    	H=String_Hash[cpg1.chrom()];
    	loci=cpg1.locus();
    	strand=cpg1.strand();
    }
    
    bool validc = false;    
    if(  ( (strand=="+" && Methy_List[H].positiveMap[methTmp] == (int)methylomes.size())  || 
    		(strand=="-" && Methy_List[H].Neg_positiveMap[methTmp] == (int)methylomes.size() ) ) && cpg1.total()>0) {
    	count_table << cpg1.chrom() << ":" << cpg1.locus() << ":" << cpg1.strand() << ":" << cpg1.context() 
                << "\t" << cpg1.meth() << "\t" << cpg1.total();
        validc=true;
    }

    meth_it++;
    while(validc && meth_it != methylomes_fstream.end()) {

	while( getline(*(*meth_it), encoding) )
	{
    		  Cpg cpg2(encoding);
                  if( (cpg2.context() != filtercontext && filtercontext != "ALL")) continue;
	    	  if( (geneid && cpg2.name()==cpg1.name()) || (!geneid && H==String_Hash[cpg2.chrom()] &&  cpg1.locus()==cpg2.locus() && cpg1.strand()==cpg2.strand() && cpg2.total()>0) )
	    	  {
      		count_table << "\t" << cpg2.meth() << "\t" << cpg2.total();
      		if(geneid) count_table << "\t" << cpg2.name();
      		break;
      	  }
  	}
      meth_it++;
    }
  }

  for(size_t ind = 0; ind < methylomes_fstream.size(); ++ind)
   	delete methylomes_fstream[ind];
}

static bool lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool lt_locus_raw_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.raw_pval < r2.raw_pval;
}

static bool ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.pos < r2.pos;
}
////BH
void adjust(vector<PvalLocus> &loci) {

      std::sort(loci.begin(), loci.end(), lt_locus_raw_pval);

      for (size_t ind = 0; ind < loci.size(); ++ind) 
      {
        const double current_score = loci[ind].raw_pval; 
        //Assign a new one.
        //const double adjust_pval 
        loci[ind].adjust_pval = loci.size()*current_score/(ind + 1); 
        //= adjust_pval;
      }

      for (vector<PvalLocus>::reverse_iterator it = loci.rbegin() + 1; it != loci.rend(); ++it) 
      {
        const PvalLocus &prev_locus = *(it - 1);
        PvalLocus &cur_locus = *(it);//cout << cur_locus.adjust_pval <<endl;
        cur_locus.adjust_pval = std::min(prev_locus.adjust_pval, cur_locus.adjust_pval);
      }

      for (vector<PvalLocus>::iterator it = loci.begin(); it != loci.end(); ++it) {
        PvalLocus &cur_locus = *(it); 
        if (cur_locus.adjust_pval > 1.0)
          cur_locus.adjust_pval = 1.0;
      }
	
      // Restore original order
      std::sort(loci.begin(), loci.end(), ls_locus_position); 
}

void fdr(vector<PvalLocus> &loci) {

      std::sort(loci.begin(), loci.end(), lt_locus_pval);

      for (size_t ind = 0; ind < loci.size(); ++ind) {
        const double current_score = loci[ind].combined_pval;

        //Assign a new one.
        const double corrected_pval = loci.size()*current_score/(ind + 1);
        loci[ind].corrected_pval = corrected_pval;
      }

      for (vector<PvalLocus>::reverse_iterator
            it = loci.rbegin() + 1; it != loci.rend(); ++it) {

        const PvalLocus &prev_locus = *(it - 1);
        PvalLocus &cur_locus = *(it);

        cur_locus.corrected_pval =
              std::min(prev_locus.corrected_pval, cur_locus.corrected_pval);
      }

      for (vector<PvalLocus>::iterator it = loci.begin();
            it != loci.end(); ++it) {
        PvalLocus &cur_locus = *(it);
	if(std::isnan(cur_locus.corrected_pval)) cur_locus.corrected_pval = 1.0;
        if (cur_locus.corrected_pval > 1.0 || cur_locus.corrected_pval < 0)
          cur_locus.corrected_pval = 1.0;
      }

      // Restore original order
      std::sort(loci.begin(), loci.end(), ls_locus_position);
}

// Splits a string using white-space characters as delimeters.
static vector<string> split(string input) {
  istringstream iss(input);
  string token;
  vector<string> tokens;

  while (iss >> token)
    tokens.push_back(token);

  return tokens;
}

// Given the maximum likelihood estimates of the full and reduced models, the
// function outputs the p-value of the log-likelihood ratio. *Note* that it is
// assumed that the reduced model has one fewer factor than the reduced model.
double loglikratio_test(double null_loglik, double full_loglik) {

  // The log-likelihood ratio statistic.
  const double log_lik_stat = -2*(null_loglik - full_loglik);

  // It is assumed that null model has one fewer factor than the full model.
  // Hence the number of degrees of freedom is 1.
  const size_t degrees_of_freedom = 1;

  // Log-likelihood ratio statistic has a chi-sqare distribution.
  double chisq_p = gsl_cdf_chisq_P(log_lik_stat, degrees_of_freedom);
  const double pval = 1.0 - chisq_p;

  return pval;
}

bool
has_low_coverage(const Regression &reg, size_t test_factor) {

  bool is_covered_in_test_factor_samples = false;
  bool is_covered_in_other_samples = false;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.design.matrix[sample][test_factor] == 1) {
      if (reg.props.total[sample] != 0)
        is_covered_in_test_factor_samples = true;
    } else {
      if (reg.props.total[sample] != 0)
        is_covered_in_other_samples = true;
    }
  }

  return !is_covered_in_test_factor_samples || !is_covered_in_other_samples;
}

bool
has_extreme_counts(const Regression &reg) {

  bool is_maximally_methylated = true;
  bool is_unmethylated = true;

  for (size_t sample = 0; sample < reg.design.sample_names.size(); ++sample) {
    if (reg.props.total[sample] != reg.props.meth[sample])
      is_maximally_methylated = false;

    if (reg.props.meth[sample] != 0)
      is_unmethylated = false;
  }

  return is_maximally_methylated || is_unmethylated;
}

string strip_methpath(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

//---------main
int main(int argc, const char **argv) {

  try {
  	  
    const string prog_name = strip_methpath(argv[0]);
    string dmc_outfile;
    string dmr_outfile;      
    double Pcutoff=0.01; double cutoff = 1; double methdiff=0.25;
    bool Auto = false;
    bool singleAuto = false;
    //unsigned length_dmr = 1000;
    unsigned methy1Start=0;
    unsigned methy1End=0;
    unsigned methy2Start=0;
    unsigned methy2End=0;
    string Genome;
    int mindmc = 4;
    int mindmcdis = 100;
    int maxdmrlen = 1000000000;
    bool gzinfile = false;
	//----help  arguments
    printf("\nBatMeth2: DMR  v1.0\n");
    const char* Help_String="Command Format :  DMR [options] -g genome.fa -o_dm <DM_result>  -1 [Sample1-methy ..] -2 [sample2-methy ..] \n"
                "\nUsage:\n"
        		"\t-o_dm        output file\n"
                "\t-g|--genome  Genome\n"
                "\t-1           sample1 methy files, sperate by space.\n"
                "\t-2           sample2 methy files, sperate by space.\n"
                "\t-o_dmr       when use auto detect by dmc.\n"
                "\t-mindmc      min dmc sites in dmr region. [default : 4]\n"
                "\t-minstep     min step in bp [default : 100]\n"
                "\t-maxdis      max length of dmr [default : 0]\n"
                "\t-pvalue      pvalue cutoff, default: 0.01\n"
                "\t-FDR         adjust pvalue cutoff default : 1.0\n"
                "\t-methdiff    the cutoff of methylation differention. default: 0.25 [CpG]\n"
                "\t-element     caculate gene or TE etc function elements.\n"
                "\t-context     Context for DM. [CG/CHG/CHH/ALL]\n"
                //"\t-f auto\n"
                "\t-L           predefinded regions or loci.\n"
                "\t-gz          gzip infile.\n"
                "\t-h|--help";
	//-----
    for(int i=1;i<argc;i++)
    {
            if(!strcmp(argv[i], "-o_dm")  )
            {
                    dmc_outfile=argv[++i];
            }
            else if(!strcmp(argv[i], "-o_dmr")  )
            {
                    dmr_outfile=argv[++i];
		    Auto=true;
                    singleAuto = true;
            }else if(!strcmp(argv[i], "-element")  )
            {
                    geneid=true;
            }else if(!strcmp(argv[i], "-context")  )
            {
                    filtercontext=argv[++i];
            }else if(!strcmp(argv[i], "-gz")  )
            {
                    gzinfile=true;
            }
            else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--genome"))
		    {
			Genome=argv[++i];
		    }
            else if(!strcmp(argv[i], "-mindmc")  )
            {
                    mindmc=atoi(argv[++i]);
            }
            else if(!strcmp(argv[i], "-minstep")  )
            {
                    mindmcdis=atoi(argv[++i]);
            }
            else if(!strcmp(argv[i], "-maxdis")  )
            {
                    maxdmrlen=atoi(argv[++i]);
            }
            else if(!strcmp(argv[i], "-FDR")  )
            {
                    cutoff=atof(argv[++i]);
            }else if(!strcmp(argv[i], "-pvalue")  )
            {
                    Pcutoff=atof(argv[++i]);
            }
            else if(!strcmp(argv[i], "-methdiff")  )
            {
                    methdiff=atof(argv[++i]);
            }
            else if( !strcmp(argv[i], "-1") )
            {
                    methy1Start=++i;
                    while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {methy1End=--i;}else {methy1End=i ;}
            }
            else if( !strcmp(argv[i], "-2") )
            {
                    methy2Start=++i;
                    while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {methy2End=--i;}else {methy2End=i ;}
            }
            else if( !strcmp(argv[i], "-L") )
            {
                    Auto=false;
		    singleAuto=false;
            }
            else
            {
                    printf("%s \n",Help_String);
                    exit(0);
            }
    }
	//
    if (argc == 1) {
	printf("%s \n",Help_String);
      return EXIT_SUCCESS;
    }
    std::string outTmp;
    if( dmc_outfile.empty() ) outTmp=dmr_outfile;
    else outTmp=dmc_outfile;
    chromLengthExact(Genome);// get chrom length
    string Genome_Len=Genome;Genome_Len+=".len";
    FILE* GenomeLen=File_Open(Genome_Len.c_str(),"r");
    
     struct Offset_Record
     {
		char Genome[40];
		unsigned Offset;
     }; 
    char Temp_OR[200]; unsigned Genome_CountX=0;
    while (fgets(Temp_OR,190,GenomeLen)!=0)//count genomes..
    {
	  Genome_CountX++;
    }
    rewind(GenomeLen);
    if(Genome_CountX==0){
        fprintf(stderr, "%s is empty! Please rerun!", Genome_Len.c_str());    
        exit(0);
    }

    Genome_Count=0;
    Offset_Record Genome_Offsets[Genome_CountX];
    while (fgets(Temp_OR,190,GenomeLen)!=0)//count genomes..
    { 
    	sscanf(Temp_OR,"%s%u",Genome_Offsets[Genome_Count].Genome,&Genome_Offsets[Genome_Count].Offset);
    	//String_Hash[Genome_Offsets[Genome_Count].Genome]=Genome_Count;
    	Genome_Count++;
    }
    
    Methy_Hash Methy_List[Genome_CountX];
    for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
    {
    		String_Hash[Genome_Offsets[i].Genome]=i;
		//Methy_List[i].positiveSample = new int[Genome_Offsets[i].Offset+1];
		//Methy_List[i].Neg_positiveSample = new int[Genome_Offsets[i].Offset+1];
		Methy_List[i].Index=i;
     }
    
    
    unsigned Sample1Size=methy1End-methy1Start+1;
    unsigned Sample2Size=methy2End-methy2Start+1;
    
   vector<string> names;
   vector<string> methylomes;
   for(unsigned i=methy1Start;i<=methy1End;i++)
   {
	names.push_back(strip_methpath( argv[i] ));
	methylomes.push_back(argv[i]);
    }
   for(unsigned i=methy2Start;i<=methy2End;i++)
   {
	names.push_back(strip_methpath( argv[i] ));
	methylomes.push_back(argv[i]);
    }
    ofstream combineTmp;
    combineTmp.open( (outTmp+"_combined.file.temp.txt").c_str() );
    ofstream combineInf;
    combineInf.open( (outTmp+"_combined.infrom.temp.txt").c_str() );
   if (!methylomes.empty()){
    if(gzinfile) merge_methylomes_gz(names, methylomes, combineTmp ,combineInf ,methy1End-methy1Start+1,Methy_List,String_Hash);
   	else merge_methylomes(names, methylomes, combineTmp ,combineInf ,methy1End-methy1Start+1,Methy_List,String_Hash);
   }

   combineInf.close();combineTmp.close();
/*
    for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
    {
		delete[] Methy_List[i].positiveSample;
		delete[] Methy_List[i].Neg_positiveSample;
     }
*/
  cout << "combined file done!\nRuning differential test ..." << endl;
  //---------------------
    // Run beta-binoimial regression or fisher test
    if(Sample1Size ==1 && Sample2Size ==1 && Auto){
	Auto=false;
	fprintf(stderr, "Without replication, Suggest use predefined region detect mode.\n");
    }

    if(Sample1Size >0 && Sample2Size >0) 
    { 
      string test_factor_name="case";
	
      ifstream design_file((outTmp+"_combined.infrom.temp.txt").c_str() );
      if (!design_file)
        throw DNAmethException("could not open file:  combined.infrom.temp.txt");
  	
	string table_filename = (outTmp+"_combined.file.temp.txt") ;
      std::ifstream table_file(table_filename.c_str());
      if (!table_file)
        throw DNAmethException("could not open file: " + table_filename);

      Regression full_regression;            // Initialize the full design
      design_file >> full_regression.design; // matrix from file.

      vector<string>::const_iterator test_factor_it =
        std::find(full_regression.design.factor_names.begin(),
                  full_regression.design.factor_names.end(), test_factor_name);

      if (test_factor_it == full_regression.design.factor_names.end())
        throw DNAmethException("Error: " + test_factor_name +
                                " is not a part of the design specification.");

      size_t test_factor = test_factor_it -
                                  full_regression.design.factor_names.begin();

      Regression null_regression;
      null_regression.design = full_regression.design;
      remove_factor(null_regression.design, test_factor);

      string sample_names_encoding;
      getline(table_file, sample_names_encoding);
      if (full_regression.design.sample_names != split(sample_names_encoding))
        throw DNAmethException(sample_names_encoding + " does not match factor "
                                "names or their order in the design matrix. "
                                "Please verify that the design matrix and the "
                                "proportion table are correctly formatted.");

      // Performing the log-likelihood ratio test on proportions from each row of the proportion table.
      static vector<PvalLocus> pvals;
      string bin_spec = "1:200:1";
      BinForDistance bin_for_dist(bin_spec);
      size_t chrom_offset = 0;std::string prev_chrom;
      while (table_file >> full_regression.props) {

        if (full_regression.design.sample_names.size() !=
            full_regression.props.total.size())
              throw DNAmethException("There is a row with"
                                      "incorrect number of proportions.");

        size_t coverage_factor = 0, coverage_rest = 0,
               meth_factor = 0, meth_rest = 0;

        for(size_t s = 0; s < full_regression.design.sample_names.size(); ++s) {
          if(full_regression.design.matrix[s][test_factor] == 0) { //!=
          	  meth_factor += full_regression.props.meth[s];
            coverage_factor += full_regression.props.total[s];
          } else {
          	  meth_rest += full_regression.props.meth[s];
            coverage_rest += full_regression.props.total[s];
          }
        }

        // Do not perform the test if there's no coverage in either all case or
        // all control samples. Also do not test if the site is completely
        // methylated or completely unmethylated across all samples.
        double pval;
        if (has_low_coverage(full_regression, test_factor)) {
          pval= -1;
        }
        else if (has_extreme_counts(full_regression)) {
          pval= -1;
        }
        else {
        	if(Sample1Size >1 && Sample2Size >1) 
        	{
	          fit(full_regression);
	          null_regression.props = full_regression.props;
	          fit(null_regression);
	          pval = loglikratio_test(null_regression.max_loglik,
	                                         full_regression.max_loglik);

	          // If error occured in the fitting algorithm (i.e. p-val is nan or -nan).
	          pval= ( (pval != pval) ? -1 : pval);
	       }else if(Sample1Size ==1 || Sample2Size ==1)
	      {
	      	  pval=fishers_exact(meth_factor/Sample1Size, coverage_factor/Sample1Size, meth_rest/Sample2Size, coverage_rest/Sample2Size) ;
	      }
	   }
	          std::string chrom=full_regression.props.chrom;
	          std::string sign = full_regression.props.strand;
	          std::string context = full_regression.props.context;
	          size_t position = full_regression.props.position;
		    std::string name=full_regression.props.name;
	          // Skip loci that do not correspond to valid p-values.
	          if(pval>1 || pval<0) pval=1;
	          if (pval >= 0 && pval <= 1) {
	            // locus is on new chrom.
	            if (!prev_chrom.empty() && prev_chrom != chrom)
	              chrom_offset += pvals.back().pos;

	            PvalLocus plocus;
	            plocus.raw_pval = pval;
	            plocus.coverage_factor=coverage_factor;
	            plocus.meth_factor=meth_factor;
	            plocus.coverage_rest=coverage_rest;
	            plocus.meth_rest=meth_rest;
	            plocus.chrom=chrom;
	            plocus.context=context;
	            plocus.position=position;
	            plocus.sign=sign;
	            plocus.name=name;
	            
	            plocus.pos = chrom_offset +
	            bin_for_dist.max_dist() + 1 + position;
                    double meth_diff=double(meth_factor)/coverage_factor - double(meth_rest)/coverage_rest ;
                    meth_diff = fabs(meth_diff);
                    if(!(meth_diff<=0.05 && pval>0.05)){
	                pvals.push_back(plocus);
			//fprintf(stderr, "%f %f\n", meth_diff, pval);
		    }
	            prev_chrom = chrom;
          }
      }
       // Combine p-values using the Z test.
    	if(pvals.size() > 0) 
	{
	      cerr << "Adjust p-values." << endl; 
	      adjust(pvals);
	      if(0 && Auto)
	      {
		      cerr << "Combining p-values." << endl;
		      combine_pvals(pvals, bin_for_dist);
		      cerr << "Running multiple test adjustment." << endl;
		      fdr(pvals);
	  	}
	}/*else if(Sample1Size ==1 && Sample2Size ==1)
	{
		fdr(pvals);
	}*/
     ofstream OutFileAdjust;
     OutFileAdjust.open(dmc_outfile.c_str());
     
     if(0 && Auto) update_pval_loci(pvals, OutFileAdjust,cutoff,methdiff);
     else 
     {
     	cout << "Printing result." <<endl;
     	Print_dm_result(pvals, OutFileAdjust,cutoff, Pcutoff, methdiff, dmr_outfile, singleAuto, mindmc, mindmcdis, maxdmrlen);
     }
//delete temp files
    remove( (outTmp+"_combined.infrom.temp.txt").c_str() );
    remove( (outTmp+"_combined.file.temp.txt").c_str() );
    }else
    {
    	printf("sample size error.");
    	exit(0);
    }
    // ------  merge DMR
    if (0 && Auto) {
     string bin_spec = "1:200:25";
     ofstream OutFiledmr;
     OutFiledmr.open(dmr_outfile.c_str());

      std::ifstream bed_file(dmc_outfile.c_str());

      if (!bed_file)
        throw "could not open file: " + dmc_outfile;
	cerr << "Definding DMR." << endl;
      merge(bed_file, OutFiledmr, Pcutoff);
      cerr << "[done]" << endl;
    } 
    /****************** END COMMAND LINE OPTIONS *****************/

  }
  catch (const DNAmethException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

