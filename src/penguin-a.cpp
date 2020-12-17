//#define NDEBUG
//Routines for pairing...  
//TODO: can we change ealignXic_MapQ to realign all? 
#define __MAIN_CODE__
#include <algorithm>

//{------------------------------  INCLUDE FILES  -------------------------------------------------/ 
#include <iostream>
#include <sstream>
#include <cstdio> 
#include "math.h"
#include <stdarg.h> 
#include <string> 
#include <getopt.h> 
#include <limits.h> 
#include <string.h> 
#include <stdlib.h> 
#include <xmmintrin.h> 
#include <emmintrin.h> 
#include <ctype.h> 
#include "bfix.h" 
#include "assert.h"
#include "ssw.h"
#include "common.h" 
#include "rqindex.h"
#include "batlib.h"
#include <map>
#include <set>
#include <queue>
#include "global.h"
#include "unistd.h"
#include "swroutines.h"
#include "print.h"
#include "filters.h"
#include <pthread.h>
#include "sched.h"
#include "fastsw.h"
#include "print.h"
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}


const unsigned MAPPED         =0;
const unsigned PE_SEQ         =1;
const unsigned PROPER_PAIR    =0x2;
const unsigned NOMAP       	 =0x4;
const unsigned MATE_UNMAPPED  =0x8;
const unsigned QUERY_MINUS    =0x10;
const unsigned MATE_MINUS     =0x20;
const unsigned FIRST_READ     =0x40;
const unsigned SECOND_READ    =0x80;
const unsigned AUX	      =0x800;

//}-----------------------------  INCLUDE FILES  -------------------------------------------------/
extern bool DO_INDEL_SMALL; //=true;
bool DO_INDEL_LARGE=true; //false;//false;
int SMALLSEED=35;
int DetectIndelLen = 100;
int Max_Hits=3000;
int time1=0,time2=0,time3=0,time4=0,time5=0;
int N1=0,N2=0,N3=0,N4=0,N5=0;
int NN1=0,NN2=0,NN3=0,NN4=0,NN_indels1=0,NN_indels2=0;
int TH1=0,TH2=0,TH2_1=0,TH3=0,TH4=0,TH5=0,TH6=0,TH7=0,TH8=0,TH9=0,TH10=0;
long File_size;
int MAXHITS_batmeth = 200;
int CLIP_SAVE_LENGTH=20;
int MIS_IN_AUX=0;
int TOP_TEN=0;
int SW_SIMILARITY_FOR_RESCUE=60;
int DISC_THRESHOLD=10;
int LENGTH_CUTOFF=0;
int BOOST=0;
int maxhits=3000;
int MAXHITS=50;//1000;
int Dummy_Int=0;
int CUT_MAX_SWALIGN=20000;
int MODE=INT_MAX;
int SEEDSIZE=INT_MAX;
int INSERT=5;
int DELETE=5;
int INSERTSIZE=800;//INT_MAX;
int STD=100;//INT_MAX;
int REAL_STD;
int SW_THRESHOLD=290;
int READS_TO_ESTIMATE=100000;
const int ScaleQ[]={0,1.5,1.75,2,3,4};
extern const Alignment Default_Alignment={0};
const int ORGSTRINGLENGTH=2200; 
bool PAIRED=FALSE;
bool Hard_Penalty=false; // false would be using --softpenalty
bool REALN=true; //false;//true;
//extern const int QUALITYCONVERSIONFACTOR=64;
//extern const int QUALITYSCALEFACTOR=33;
bool DASH_DEL=true;
bool ESTIMATE=false; //true;//
bool FASTDECODE=false;
bool DEB=false;
bool FASTSW=true;
bool DEBUG_SEGS=false;
unsigned GENOME_SIZE=0;
int QUALITYCONVERSIONFACTOR=33;
int QUALITYSCALEFACTOR=1;
BATPARAMETERS BP;
int THREAD = 8;
int Top_Penalty;
int JUMP=8;
int INDELGAP=21;//15 good for //17 for 8, 19=g00d for 9
FMFILES FMFiles,CT,GA;
bool Print_hits=false;//shi fou yi jing shu chu
typedef std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> ALIGNMENT_Q;
inline void really_free(std::map<unsigned,Alignment>& to_clear)
{
    std::map<unsigned,Alignment> v;
    v.swap(to_clear);
}
template <typename T>
inline void really_free_M(std::vector<T>& to_clear)
{
    std::vector<T> v;
    v.swap(to_clear);
}
int Nindel=0;
//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
void get_len_match(Alignment & A, int & headM, int & TailM);
bool Align_comp(Alignment & A, Alignment & B, int Aheadlen, int Ataillen, int Bheadlen, int Btaillen, int Real_len);
inline unsigned char Hash(char* S);
void Build_Names(const char* Genome_Name,FMFILES & F,BATPARAMETERS & BP);
inline void ReplaceCtoT(READ & R);
inline void ReplaceGtoA(READ & R);
unsigned ulabs(unsigned A,unsigned B);
void Remove_Dup_Final(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut, int & secondN);
int Get_ED(std::string & S);
void Map_One_SEG_Head(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE);
void Map_One_SEG_Tail(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF2,MEMX & MC2,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE);
void Map_One_Half_Head(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLH,MEMX & MCLH,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score);
void Map_One_Half_Tail(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLT,MEMX & MCLT,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score);
void Two_Side_Hit_Finding(int & align_mismatch,int & Paired_Ncutoff,char source1,char source2,READ & RawR,READ & RawM,int Read_Length,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,READ & M,BATREAD & B,unsigned & Conversion_Factor,MEMX & mF,MEMX & mC,MEMX & MF,MEMX & MC,MEMX & MF2,MEMX & MC2,MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Main,LEN & L_Half,LEN & L_Third,unsigned & Actual_Tag,FILE*  Single_File,FILE* Mishit_File,Align_Hit & Align_Hits,Align_Hit & Align_Hits_P,Pair* & Pairs,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut);
void Map_One_Quart_Head(int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Third,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int SEG_SIZE,int SHIFT_SEG_SIZE,Pair* & Pairs);
void Map_One_Quart_Tail(int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF2,MEMX & MC2,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Third,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int SEG_SIZE,int SHIFT_SEG_SIZE,Pair* & Pairs);

void MEM_STRUC(MEMX & MF, MEMX & MC,MEMX & MF2, MEMX & MC2,MEMX & MCLH,MEMX & MFLH, MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,MEMLOOK & MLook,LEN & L,MEMX & mF,MEMX & mC);
void FileInfo(char *PATTERNFILE,char *HITSFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID);
void Mapping_Batmis(int & align_mismatch,int & Paired_Score,char source1,char source2,READ & R,READ & M,BATREAD & B,MEMX & MF,MEMX & MC,MEMX & MF2,MEMX & MC2,MEMX & mF,MEMX & mC,LEN & L_batmeth,BWT* fwfmi,BWT* revfmi,READ & RawR,READ & RawM,Align_Hit & Align_Hits,Align_Hit & Align_Hits_P,OUTPUT & O1,OUTPUT & O2,unsigned char* Original_Text,FILE*  Single_File,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,
		MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Main,LEN & L_Half,LEN & L_Third,unsigned & Actual_Tag,FILE* Mishit_File,Pair* & Pairs,unsigned & Conversion_Factor,unsigned Entries,RQINDEX & RQHALF,RQINDEX & RQ);
bool Output_Pair(Alignment A1,Alignment A1_P,Alignment B1,Alignment B1_P,int Read_Length,int Read_Length2);
void Extend_Right(char* Temp_Current_Tag_raw,READ & RawR,RQINDEX & RQHALF,unsigned char* Original_Text,int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter );
void Extend_Left(char* Temp_Current_Tag_raw,READ & RawR,RQINDEX & RQHALF,unsigned char* Original_Text,int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter);
int Do_Indel(char source,READ & RawR,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT *fwfmi,BWT *revfmi,MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,int StringLength,Pair* & Pairs,FILE* & Single_File,READ & R,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,const Hit_Info & H);
void Show_Progress(unsigned Percentage);
void Read_INI(char* Config_File,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP);
void Get_Bases (unsigned char* Original_Text,unsigned Location,int StringLength,char* Org_String);
void Pair_Reads(BWT* fwfmi,BWT* revfmi,RQINDEX & R,FILE* Data_File, In_File IN,unsigned MAXCOUNT,BATPARAMETERS BP,Pair* Pairs,unsigned Entries);
void Get_Best_Alignment_Pair(READ & RawR,char source,unsigned char* Original_Text,Alignment & A,Alignment & B,READ & R,const int StringLength,BATREAD & Read,Hit_Info & H,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,int & Clip_H,int & Clip_T,char* CIG,bool PRINT,bool DUMMY_FORCED);
void Parse_Command_line(int argc, char* argv[],char* & GENOME,unsigned & MAXCOUNT,FMFILES & FMFiles,BATPARAMETERS & BP);
int Head_Tail(BWT* revfmi,SARange* Head_Hits,SARange* Tail_Hits,int Insert,int MAXCOUNT,char & In_Large,RQINDEX & R,unsigned Entries,Pair* Pairs,int & Pairs_Index,int & HITS,int & Err,unsigned Conversion_Factor);
void *Map_And_Pair_Solexa(void *T);
void Verbose(char *BWTFILE,char *OCCFILE,char *REVBWTINDEX,char *REVOCCFILE,char *PATTERNFILE,char *HITSFILE,char* LOCATIONFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID);
void Init(In_File & IN,FMFILES F,RQINDEX R,BATPARAMETERS & BP,char Solid,char Npolicy,LEN & L);
void  Paired_Extension(char* Fwd_Read_raw,char *Revcomp_Read_raw,READ & RawR,unsigned char* Original_Text,unsigned Entries,BWT *revfmi,int Last_MisT,int Last_MisH,char* Fwd_Read,char *Revcomp_Read, RQINDEX & RQ,Pair* & Pairs,SARange* & MFH_Hit_Array,SARange* & MFT_Hit_Array,SARange* & MCH_Hit_Array,SARange* & MCT_Hit_Array,int StringLength,READ & R,int & Err,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Tot_SW_Scans);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
bool Report_SW_Hits(READ & RawR,char source,unsigned char* Original_Text,const int Err,READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,Hit_Info & Mismatch_Hit,int Quality_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,bool PRINT,bool DUMMY_FORCED=false);
bool Unique_Hits(int Plus_Hits,int Minus_Hits,SARange & P,SARange & M);
bool Check_Subopt(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,float & Top_BQScore,float & Sub_Score,float & Sub_BQScore, int & Quality_Score);
Alignment Realign(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength,READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments);
Alignment Realign(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength,READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments);
Alignment RealignX(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength,READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,char* Cigar,int & Clip_H,int & Clip_T,int & Filter=Dummy_Int,bool Do_Filter=false);
bool Remove_Duplicates_And_Filter_Top_Hits(Hit_Info *Hit_List,int & Hit_Count,int StringLength);
float Calc_Top_Score(MEMX & MF,MEMX & MC,float & Top_BQ,int Top_Mis,int StringLength,int Plus_Hits,int Minus_Hits,READ & R);
unsigned SA2Loc(BWT *revfmi,SARange S,int Pos,unsigned Conversion_Factor);
void Mode_Parameters(BATPARAMETERS BP);
void Set_Force_Indel(bool & Force_Indel,int Last_Mis,Hit_Info & H,int Avg_Q);
int Calculate_Average_Quality(READ & R);
void Set_Affinity();
Alignment RealignFast(READ & RawR,char source,unsigned char* Original_Text,Hit_Info &  H,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter);
void FreeQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t );
void Get_top1M_Alignments(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t_n,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t, int & cutoff ) ;
void FreeQ_Pair(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair>  & t );
void FreeQ_Hits(Align_Hit & Align_Hits);
bool Report_Single(READ & RawR,unsigned char* Original_Text,READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,bool & Print_Status,int Clip_H,int Clip_T,Alignment & A);
void Get_Basic_MapQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Alignment & C1, int & MapQ);
void Adjust_Alignments(READ & RawR,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp);
void Pop_And_Realign(READ & RawR,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,int Offset,READ & RTemp, BATREAD & BTemp);
bool SW_List(READ & RawR,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp,bool & List_Exceeded,int & List_Size);
bool Correct_Orientation(Alignment A,Alignment A_P,int Extra_Bit=0);
bool Find_Paired(int & paired_score,bool & Unique,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & B,std::map<unsigned,Alignment> & D,std::map<unsigned,Alignment> & D_P,int Extra_Bit=0,int SW_Compare=0);
bool Align_Difference(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,unsigned U);
void Mate_Rescue_New(int mismatch,int Indel,bool & Do_Rescue,bool & realign,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,READ & RawR,READ & RawM,unsigned char* Original_Text,char source1,char source2,READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,int Read_Length,int Read_Length2,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments_P,Hit_Info & H1,Hit_Info & H1_P,FILE* Single_File,int Quality_Score1,int Quality_Score1_P,Alignment & A1,int MapQ2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi);
void Remove_Dup_Top(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,int Gap);
void Fix_Offset(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off);
void Fix_Offset_Head(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off);
void Fix_Offset_Tail(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off);
void Print_Pair(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,Hit_Info &  HitsH,Hit_Info &  HitsT,char Sign1,char Sign2,READ & RawR,READ & RawM,unsigned char* Original_Text,char source1,char source2,FILE* Single_File,Final_Hit & H,Final_Hit & T,READ & R1, READ & R2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi);
int Get_Skip(std::string & S);
bool Check_Proper_Pair(int S1,int S2,unsigned Loc1, unsigned Loc2,int Extra_Bit);
void Estimate_Insert(int & INSERTSIZE,int & STD);
void Rescue_Clip(std::string S,int & P_Clip,int & S_Clip);
void Print_Aux_Hit(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,char Sign,READ & RawR,char source,FILE* Single_File,Final_Hit & H,Final_Hit & Aux,READ & R,int Clip1,int Clip2,MEMX & MF);
pthread_mutex_t Lock_Estimate;
std::vector<int> Estimate;
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
std::string GetFilePosfix(const char* path)
{
//    std::cout << strlen(path) << std::endl;
    std::string str(path);
    str=str.substr(str.length()-2,2);
    std::transform(str.begin(),str.end(),str.begin(),::tolower);   
    return str;    
}

void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
  std::string::size_type pos1, pos2;
  pos2 = s.find(c);
  pos1 = 0;
  while(std::string::npos != pos2)
  {
    v.push_back(s.substr(pos1, pos2-pos1));

    pos1 = pos2 + c.size();
    pos2 = s.find(c, pos1);
  }
  if(pos1 != s.length())
    v.push_back(s.substr(pos1));
}

#undef DEBUG
//bool TESTMODE=false;
const bool EXACT=false;
const int QUAL_UNMAPPED=20000;
int SW_STRING_BUFFER=1400;
const char* Code_To_CharCAPS="ACGT";
const int DISK_BUFFER_SIZE=5000;
const int MAXALN=30000;
const int MAXCOUNT=200;
const int MIN_SCORE_DIFF=10;
bool PRINT_ALL_HITS=false;
int INDEX_RESOLUTION=30000;
int Inter_MM = 5;
int Max_MM_GAP = 2;
int Max_MM_GAP_Adjust = 2;
int Max_MM = 5;
int MAX_SW_HITS = 100;
int MAX_SW_sav;
bool USE_MULTI_OUT=false ;//true;//true;
bool PRINT_HEADER=true; //false;
bool MAIN_HEADER_ONLY=false;;
bool STACK_LOWQ=false;
const int DEFAULTFASTAQUAL=35;//40;
LEN L_Main;
In_File IN;
//INFILE Head_File;
INFILE Tail_File;
time_t Start_Time,End_Time;
FILE* Main_Out;
int gap_openP,gap_extensionP;
int MAX_PER_LIST=60;
int SCOREGAP=20;
std::string RGID;
extern char Char_To_R[255];
int match_SC=1;
extern bool non_directional;
MEMLOOK MLookGA,MLookCT;

int Genome_CountX=0;
Offset_Record* Genome_Offsets;
std::map <std::string,int> String_Hash;
bool gzinfile = false;
std::vector <std::string> filelist1;
std::vector <std::string> filelist2;
pthread_mutex_t flists;
std::vector <INFILE> Ffilelist1;
std::vector <INFILE> Ffilelist2;
std::vector <INFILE>::size_type fi = 0;
int processed = -1;
//{---------------------------- GLOBAL VARIABLES -------------------------------------------------


int main(int argc, char* argv[])
{
	FILE* Data_File;
	FILE* Mishit_File;
	char *CT_GENOME,*GA_GENOME,*GENOME;//moxian
	MMPool *mmPool;
	unsigned MAXCOUNT=1;
	BP.PAIRING_MODE=NORMALFILEMODE;BP.NTHREADS=0;BP.FORCELENGTH=0;BP.ROLLOVER=FALSE;BP.SCANBOTH=FALSE;BP.ONEFMINDEX =FALSE; BP.MAXHITS=1;BP.CMD_Buffer=new char [5000];

	gap_openP=40,gap_extensionP=6;
	//Read_INI(NULL,MAXCOUNT,FMFiles,BP);
	//Read_INI(NULL,MAXCOUNT,CT,BP);
	Read_INI(NULL,MAXCOUNT,GA,BP);
	Parse_Command_line(argc,argv,GENOME,MAXCOUNT,GA,BP);
	//------------------------------CT GA genome--------------
	std::string CTS=GENOME,GAS=GENOME;
	if (GENOME)
	{
		CTS+="-CtoT";GAS+="-GtoA";//  -CtoT -CtoT-GtoA
		CT_GENOME=(char*)CTS.c_str();GA_GENOME=(char*)GAS.c_str();
	}

	Build_Names(GENOME,FMFiles,BP);
	Build_Names(GA_GENOME,GA,BP);// FMFiles--GA
	Build_Names(CT_GENOME,CT,BP);

        //-----------------------------batmeth1
        std::string L=GENOME;L+=".ann.location";

        FILE* Location_File=File_Open(L.c_str(),"r");

        char Genome[40];
        while (fgets(Genome,39,Location_File)!=0)//count genomes..
        {
                fgets(Genome,39,Location_File);
                Genome_CountX++;
        }
        rewind(Location_File);
        Genome_Offsets = new Offset_Record[Genome_CountX+1];

        int Genome_Count=0;
        while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<=Genome_CountX)
        {
                Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
                fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
                for(int i=0;i<40;i++)
                {
                        if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
                        {
                                Genome_Offsets[Genome_Count].Genome[i]=0;
                                break;
                        }
                }
                String_Hash[Genome_Offsets[Genome_Count].Genome]=Genome_Count;
                Genome_Offsets[Genome_Count].index=Genome_Count;
                Genome_Count++;
        }
        Genome_Count--;

	//---------------------------------------------------------------	
	init_SSW();Build_Pow10();
	init_SSW_Clip(match, mismatch, gap_open, gap_extension);
	//init_SSW_Clip(match_SC /*match*/,3 /*mismatch*/,11 /*gap_open*/,4 /*gap_extension*/);
	if(CONFIG_FILE){ Read_INI(CONFIG_FILE,MAXCOUNT,GA,BP); Read_INI(CONFIG_FILE,MAXCOUNT,CT,BP); }
	if (MISC_VERB) fprintf(stderr,"BatMeth2 v2.00\n");
	if (BP.MAX_MISMATCHES != INT_MAX)
	{
		if (BP.MAX_MISMATCHES <=5) Inter_MM=BP.MAX_MISMATCHES;
		else fprintf(stderr,"Error tolerence too high ...\n");
	}
	

	Match_FA_Score= -10*log10(Pr(DEFAULTFASTAQUAL));
	Mis_FA_Score= -10*log10((1-Pr(DEFAULTFASTAQUAL))/3);
	SplitString(BP.PATTERNFILE, filelist1, ",");
        if(PAIRED) SplitString(BP.PATTERNFILE1, filelist2, ",");
        std::string filefm = GetFilePosfix(filelist1[0].c_str());
        if(filefm=="gz") gzinfile=true;

	for(std::vector <std::string>::size_type i = 0;i<filelist1.size();i++){
		INFILE Head_File;
	        if(gzinfile){
        	        Head_File.gzfp = gzopen(filelist1[i].c_str(), "rb");
	                if(PAIRED) Tail_File.gzfp = gzopen(filelist2[i].c_str(), "rb");
                	Head_File.filegz = Tail_File.filegz = '1';
        	}else{
	                Head_File.Input_File=File_Open(filelist1[i].c_str(),"r");
                	if(PAIRED) Tail_File.Input_File=File_Open(filelist2[i].c_str(),"r");
        	        Head_File.filegz = Tail_File.filegz = '0';
	        }
		Analyze_File(Head_File,L_Main);
		Ffilelist1.push_back(Head_File);
	}
	IN.HEAD_LENGTH=IN.TAIL_LENGTH=IN.STRINGLENGTH=L_Main.STRINGLENGTH;

	//Load_Range_Index(RQ,L_Main.STRINGLENGTH/3,FMFiles,Entries);
	Load_Range_Index(RQ_GA,L_Main.STRINGLENGTH/3,GA,Entries_GA);//
	Load_Range_Index(RQ_CT,L_Main.STRINGLENGTH/3,CT,Entries_CT);
	if(FASTDECODE) {Load_Range_Index(RQHALF_GA,L_Main.STRINGLENGTH/2,GA,Entries_Half_GA);Load_Range_Index(RQHALF_CT,L_Main.STRINGLENGTH/2,CT,Entries_Half_CT);}
	
	//Init(IN,FMFiles,RQ,BP,Head_File.SOLID,0,L_Main);
	Init(IN,GA,RQ_GA,BP,Ffilelist1[0].SOLID,0,L_Main);
	fprintf(stderr, "==============================================================\nLoading genome: %s ...\n",GENOME);
	//Init(IN,CT,RQ_CT,BP,Head_File.SOLID,0,L_Main); 
	//---------------------------------------------------------------------------------------------------------------------------------
	//--remove
	FILE* Original_File=File_Open(FMFiles.BINFILE,"rb");
	Original_Text_Ori=(unsigned char*) malloc(Get_File_Size(Original_File));
	if(!fread(Original_Text_Ori,Get_File_Size(Original_File),1,Original_File))fprintf (stderr,"Init(): Error loading genome...\n");


	FILE* Original_File_GA=File_Open(GA.BINFILE,"rb");
	Original_Text_GA=(unsigned char*) malloc(Get_File_Size(Original_File_GA));
	if(!fread(Original_Text_GA,Get_File_Size(Original_File_GA),1,Original_File_GA))fprintf (stderr,"Init(): Error loading genome...\n");
	File_size=sizeof(unsigned)*Get_File_Size(Original_File_GA);
	
	FILE* Original_File_CT=File_Open(CT.BINFILE,"rb");
	Original_Text_CT=(unsigned char*) malloc(Get_File_Size(Original_File_CT));
	if(!fread(Original_Text_CT,Get_File_Size(Original_File_CT),1,Original_File_CT))fprintf (stderr,"Init(): Error loading genome...\n");
	//-----------------------------------------------------------------------------------------------------------------------------------
	Mode_Parameters(BP);

	unsigned Location_Array[80];

	Load_FM_and_Sample(fwfmiGA,revfmiGA,mmPool,GA); 
	Load_FM_and_Sample(fwfmiCT,revfmiCT,mmPool,CT);


	//---------------------Load Genome-----------------------------------
	//-------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------	
	GENOME_SIZE=revfmiGA->textLength;assert(GENOME_SIZE==fwfmiGA->textLength);
	Load_Location(GA.LOCATIONFILE,Annotations,S,E,Location_Array);
	if (NFILE) Load_LocationN(GA.NLOCATIONFILE,Annotations,S,E,Location_Array);
	if(!USE_MULTI_OUT)
	{
		Main_Out=File_Open(GA.OUTPUTFILE,"w");
		if(PRINT_HEADER)
		{
			std::map <unsigned, Ann_Info> ::iterator S,E;
			S=Annotations.begin();E=Annotations.end();
			unsigned CSize=0;
			while (S!=E)
			{
				Ann_Info T=S->second;
				if(T.Size > CSize) CSize=T.Size;
				fprintf(Main_Out,"@SQ\tSN:%s\tLN:%d\n",T.Name,T.Size);
				S++;
			}
			if(NFILE && CSize) fprintf(Main_Out,"@SQ\tSN:%s\tLN:%d\n","NOISE",CSize);
			char Current_Dir[1000];
			if (!getcwd(Current_Dir,990))
			{
				sprintf (Current_Dir,"%d",rand());
			}
			std::ostringstream ostr;
			ostr << (rand()%(INT_MAX));
			RGID=ostr.str();
			if(PAIRED)
            	fprintf(Main_Out,"@RG\tID:%s\tSM:%s\tLB:%s,%s\n",RGID.c_str(),Current_Dir,BP.PATTERNFILE,BP.PATTERNFILE1);
            else fprintf(Main_Out,"@RG\tID:%s\tSM:%s\tLB:%s\n",RGID.c_str(),Current_Dir,BP.PATTERNFILE);
			fprintf(Main_Out,"@PG\tID:PEnGuin\tCL:%s",BP.CMD_Buffer);
		}
	}
		// reverse
	Char_To_R['A']='T';Char_To_R['a']='t';
	Char_To_R['C']='G';Char_To_R['c']='g';
	Char_To_R['G']='C';Char_To_R['g']='c';
	Char_To_R['T']='A';Char_To_R['t']='a';
	Char_To_R['N']='N';Char_To_R['n']='n';
	// Initialise indexes 
	//
	MLookCT.Lookupsize=MLookGA.Lookupsize=3;//Get_Lookup_Size(BP.MAX_MISMATCHES,L.STRINGLENGTH);
		Build_Tables(fwfmiGA,revfmiGA,MLookGA);
		Build_Tables(fwfmiCT,revfmiCT,MLookCT);
//********************************************************************************************************
	Thread_Arg T;
	if (THREAD)
	{
		//fprintf(stderr,"Estimating insert size..\n");
		READS_TO_ESTIMATE=0;//READS_TO_ESTIMATE/THREAD;
		//Launch_Threads(THREAD, Map_And_Pair_Solexa,T);
		//Estimate_Insert(INSERTSIZE,STD);
		Total_Reads = 0; Total_Mapped = 0; Total_Paired = 0;
		ESTIMATE=false;
		fi=0;
		time(&Start_Time);
		Launch_Threads(THREAD, Map_And_Pair_Solexa,T);
	}
	else
	{
		Set_Affinity();
		ESTIMATE=false;
		time(&Start_Time);
		Map_And_Pair_Solexa(NULL);
	}

//********************************************************************************************************
    if(gzinfile) fprintf_time(stderr, "Processed %d reads.\n", Total_Reads);
    else Done_Progress();
	//if (PROGRESSBAR) fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
	//if(PROGRESSBAR) fprintf(stderr, "\r[+++++++++++++++++++++++++++++++++++++++100%++++++++++++++++++++++++++++++++++++++++++++++++++]\n");
	if (MISC_VERB)
	{
		fprintf(stderr,"%d / %d Reads \n",Total_Mapped,Total_Reads);
		//alignresults.txt
		fprintf(Log_SFile,"Total_Reads\t%u\n",Total_Reads);
		fprintf(Log_SFile,"Total_Mapped\t%u\n",Total_Mapped);
		fprintf(Log_SFile,"MappedQ20\t%u\n",passed20);
		//printf("%d Large reads....\n",Large);
		time(&End_Time);fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	//if (LOG_SUCCESS_FILE) fprintf(Log_SFile,"DONE\n");
}

//------------------------------- Print /Verify the reads ----------------------------------------------------
void *Map_And_Pair_Solexa(void *T)
{
	Thread_Arg *TA;
	TA=(Thread_Arg*) T;
	int Thread_ID;
	if(T)
	{
		Thread_ID=TA->ThreadID;
	}

	Align_Hit Align_Hits_CT;
	Align_Hit Align_Hits_CT_P;
//	Align_Hit Align_Hits_CT;
//	Align_Hit Align_Hits_CT_P;
	//--store result
	std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> Alignments_Reslut;
	//
	Pair* Pairs;
	//FILE* Output_File;
	FILE* Discordant_File;
	FILE* Single_File;
	FILE* Unmapped_File;
	FILE* Multi_File;
	FILE* Mishit_File;
	MMPool *mmPool;
	int Half=0;
	LEN L=L_Main,L_Half,L_Third;L.IGNOREHEAD=0;
	OUTPUT O;
	HEADER Header;
	//time_t Start_Time,End_Time;
	//time(&Start_Time);
	int MF1_Top_End,MC1_Top_End;
	char* Disc_Hit_Buffer;
#ifndef NDEBUG
	//DISC_HIT_BUFFER_ORG=Disc_Hit_Buffer_Org;
#endif
//{--------------------------- INIT STUF ---------------------------------------
	LEN L_batmeth=L;
	L_Half=L_Third=L;
	LEN mL=L;
	if (BP.FORCELENGTH) 
	{
		Split_Read(BP.FORCELENGTH,L);SEEDSIZE=BP.FORCELENGTH;
	}
	else 
	{
		Split_Read(L.STRINGLENGTH_ORG,L);SEEDSIZE=L.STRINGLENGTH_ORG;
	}
	if (!(Pairs=(Pair*)malloc(sizeof(Pair)*30000))) {fprintf(stderr,"Init():malloc error...\n");exit(-1);}
	Split_Read(Half,L_Half);
	Split_Read(IN.STRINGLENGTH/3,L_Third);
	Split_Read(20,mL);

// Misc stuff 
	if(THREAD==0 || THREAD==1) {Thread_ID=1;}
	if(Thread_ID==1) 
	{
		//Verbose(CT.BWTFILE,CT.OCCFILE,CT.REVBWTINDEX,CT.REVOCCFILE,BP.PATTERNFILE,CT.OUTPUTFILE,CT.LOCATIONFILE,Inter_MM,BP.Patternfile_Count,BP.PATTERNFILE1,Head_File.FILETYPE,L,BP.FORCESOLID);
		//Verbose(GA.BWTFILE,GA.OCCFILE,GA.REVBWTINDEX,GA.REVOCCFILE,BP.PATTERNFILE,GA.OUTPUTFILE,GA.LOCATIONFILE,Inter_MM,BP.Patternfile_Count,BP.PATTERNFILE1,Head_File.FILETYPE,L,BP.FORCESOLID);
		FileInfo(BP.PATTERNFILE,GA.OUTPUTFILE,Inter_MM,BP.Patternfile_Count,BP.PATTERNFILE1,Ffilelist1[0].FILETYPE,L,BP.FORCESOLID);
	}
	if(!USE_MULTI_OUT)
	{
		Single_File=Main_Out;
	}
	else
	{
		if(THREAD==0 || THREAD==1) {Single_File=File_Open(GA.OUTPUTFILE,"w");Thread_ID=1;}
		else
		{
			std::string Temp=GA.OUTPUTFILE;
			char Temp_Char[20];sprintf(Temp_Char,".%d",Thread_ID);
			Temp+=Temp_Char;
			Single_File=File_Open(Temp.c_str(),"w");
		}
	}

	if(PRINT_MISHIT) Mishit_File=File_Open(BP.MISFILE1,"w");
	if (!ESTIMATE && PRINT_DICTIONARY && USE_MULTI_OUT && PRINT_HEADER)
	{
		std::map <unsigned, Ann_Info> ::iterator S,E;
		S=Annotations.begin();E=Annotations.end();
		unsigned CSize=0;
		while (S!=E)
		{
			Ann_Info T=S->second;
			if(T.Size > CSize) CSize=T.Size;
			fprintf(Single_File,"@SQ\tSN:%s\tLN:%d\n",T.Name,T.Size);
			S++;
		}
		if(NFILE && CSize) fprintf(Single_File,"@SQ\tSN:%s\tLN:%d\n","NOISE",CSize);
		char Current_Dir[1000];
		if (!getcwd(Current_Dir,990))
		{
			sprintf (Current_Dir,"%d",rand());
		}
		std::ostringstream ostr;
		ostr << (rand()%(INT_MAX));
		RGID=ostr.str();
		fprintf(Single_File,"@RG\tID:%s\tSM:%s\tLB:%s\n",RGID.c_str(),Current_Dir,BP.PATTERNFILE);
		fprintf(Single_File,"@PG\tID:PEnGuin\tCL:%s",BP.CMD_Buffer);
	}
	if (WRITE_DISCORDANT) Discordant_File=fopen(GA.DISCORDANTFILE,"w");else Discordant_File=Single_File;
	//if (WRITE_SINGLE) Single_File=fopen(FMFiles.SINGLEFILE,"w");else 
	if (WRITE_UNMAPPED) Unmapped_File=fopen(BP.UNMAPPED_FILE,"w");else Unmapped_File=Single_File;
	if (WRITE_MULTI) Multi_File=fopen(BP.MULTI_FILE,"w");else Multi_File=Single_File;
	WRITE_DISCORDANT=WRITE_MULTI=WRITE_UNMAPPED=WRITE_SINGLE=TRUE;

	//FILE* Log_File=File_Open(LOGFILE,"w");GLog_File=Log_File;
	unsigned Mapped=0,Actual_Tag=0,Large=0; //,Total_Reads=0;
	int ntemp=0;Nindel=0;
	int LOOKUPSIZE; //MAX_MISMATCHES=BP.MAX_MISMATCHES;
	char ONEFMINDEX=BP.ONEFMINDEX;
	READ R;BATREAD B;BATREAD B_batmeth;
	READ revR,R_CT,R_GA;//+ strand reads..
	READ M;
	READ revM,M_CT,M_GA;//+ strand reads..
	GUESS G;GUESS G1;
	//------GA--------
	MEMX MF_GA,MC_GA,MF_GA2,MC_GA2,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA;//MEMX MF1,MC1;
	MEMX MFH_GA,MCH_GA,MFT_GA,MCT_GA;//One third pref/suf..
	MEMX mF_GA,mC_GA;
	//----CT----
	MEMX MF_CT,MC_CT,MF_CT2,MC_CT2,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT;//moxian
	MEMX MFH_CT,MCH_CT,MFT_CT,MCT_CT;
	MEMX mF_CT,mC_CT;
	
// calculate read portions 
	LOOKUPSIZE=3;//Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	//--------------------------------------GA----------------------------------
	MF_GA.Lookupsize=LOOKUPSIZE;MC_GA.Lookupsize=LOOKUPSIZE;MFLH_GA.Lookupsize=LOOKUPSIZE;MCLH_GA.Lookupsize=LOOKUPSIZE;
	MFLT_GA.Lookupsize=LOOKUPSIZE;MCLT_GA.Lookupsize=LOOKUPSIZE;MF_GA2.Lookupsize=LOOKUPSIZE;MC_GA2.Lookupsize=LOOKUPSIZE;
	MFH_GA.Lookupsize=LOOKUPSIZE;MCH_GA.Lookupsize=LOOKUPSIZE;MFT_GA.Lookupsize=LOOKUPSIZE;MCT_GA.Lookupsize=LOOKUPSIZE;
	MF_GA.L=MC_GA.L=MF_GA2.L=MC_GA2.L=MFLH_GA.L=MCLH_GA.L=MFLT_GA.L=MCLT_GA.L=MFH_GA.L=MCH_GA.L=MFT_GA.L=MCT_GA.L=L;
	mF_GA.L=mC_GA.L=mL;
	mF_GA.Lookupsize=LOOKUPSIZE;mC_GA.Lookupsize=LOOKUPSIZE;
	MF_GA.batmeth1=MC_GA.batmeth1=MF_GA2.batmeth1=MC_GA2.batmeth1=MFLH_GA.batmeth1=MCLH_GA.batmeth1=MFLT_GA.batmeth1=MCLT_GA.batmeth1=MFH_GA.batmeth1=MCH_GA.batmeth1=MFT_GA.batmeth1=MCT_GA.batmeth1=false;
	mF_GA.batmeth1=mC_GA.batmeth1=false;
	//--------------------------------------CT-------------------------------------------
	MF_CT.Lookupsize=LOOKUPSIZE;MC_CT.Lookupsize=LOOKUPSIZE;MFLH_CT.Lookupsize=LOOKUPSIZE;MCLH_CT.Lookupsize=LOOKUPSIZE;
	MFLT_CT.Lookupsize=LOOKUPSIZE;MCLT_CT.Lookupsize=LOOKUPSIZE;MF_CT2.Lookupsize=LOOKUPSIZE;MC_CT2.Lookupsize=LOOKUPSIZE;
	MFH_CT.Lookupsize=LOOKUPSIZE;MCH_CT.Lookupsize=LOOKUPSIZE;MFT_CT.Lookupsize=LOOKUPSIZE;MCT_CT.Lookupsize=LOOKUPSIZE;
	MF_CT.L=MC_CT.L=MF_CT2.L=MC_CT2.L=MFLH_CT.L=MCLH_CT.L=MFLT_CT.L=MCLT_CT.L=MFH_CT.L=MCH_CT.L=MFT_CT.L=MCT_CT.L=L;
	mF_CT.L=mC_CT.L=mL;mF_CT.Lookupsize=LOOKUPSIZE;mC_CT.Lookupsize=LOOKUPSIZE;
	MF_CT.batmeth1=MC_CT.batmeth1=MF_CT2.batmeth1=MC_CT2.batmeth1=MFLH_CT.batmeth1=MCLH_CT.batmeth1=MFLT_CT.batmeth1=MCLT_CT.batmeth1=MFH_CT.batmeth1=MCH_CT.batmeth1=MFT_CT.batmeth1=MCT_CT.batmeth1=false;
        mF_CT.batmeth1=mC_CT.batmeth1=false;
// initialise memory structures 
	MEM_STRUC(MF_CT,MC_CT,MF_CT2,MC_CT2,MCLH_CT,MFLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT,MLookCT,L,mF_CT,mC_CT);
	MEM_STRUC(MF_GA,MC_GA,MF_GA2,MC_GA2,MCLH_GA,MFLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,MLookGA,L,mF_GA,mC_GA);

//------------------batmeth1
/*      MEMX MF_CT_bat,MC_CT_bat;
        MEMX MF_GA_bat,MC_GA_bat;
        MEMX mF_CT_bat,mC_CT_bat;
        MEMX mF_GA_bat,mC_GA_bat;
*/
	L_batmeth.IGNOREHEAD=0;L_batmeth.STRINGLENGTH=70;
        L_batmeth.STRINGLENGTH_ORG=L.STRINGLENGTH;
        //Initial_Length(L_batmeth);
        Split_Read(L_batmeth.STRINGLENGTH_ORG,L_batmeth);
        B_batmeth.IGNOREHEAD=L_batmeth.IGNOREHEAD;B_batmeth.StringLength=L_batmeth.STRINGLENGTH;

//------------------batmeth
        OUTPUT O_CT,O_GA,O_RevCT,O_RevGA;

        O_CT.SAM=0;
        O_CT.PLUSSTRAND=0;
        O_CT.MaxHits=MAXHITS_batmeth;
        O_GA.MaxHits=MAXHITS_batmeth;
        O_CT.Offset=0;
        O_CT.Genome_Count=Genome_CountX;
        O_CT.FILETYPE=FQ;
        O_CT.Length_Array[0]=O_CT.Length_Array[1]=O_CT.Length_Array[2]=L_batmeth.STRINGLENGTH;
        O_RevCT=O_CT;
        O_GA=O_CT;O_RevGA=O_CT;
	O_GA.Offset=0;
        B_batmeth.IGNOREHEAD=L_batmeth.IGNOREHEAD;B_batmeth.StringLength=L_batmeth.STRINGLENGTH;
//-------------------------

	int Progress=0;unsigned Number_of_Tags=1000;
	int Bindel=0;
	if (!gzinfile && PROGRESSBAR && !ESTIMATE) Init_Progress();
	int NCount_For_SW=L.STRINGLENGTH/4;
	unsigned Conversion_Factor;

	bwase_initialize(); 
	INPUT_FILE_TYPE=Ffilelist1[0].FILETYPE;
	MAX_SW_sav=MAX_SW;
	int Read_Length=Ffilelist1[0].TAG_COPY_LEN;
	SW_THRESHOLD=80*Read_Length*match/100;
//}--------------------------- INIT STUF ---------------------------------------
    for(; fi < filelist1.size(); ){
		if(!ESTIMATE && Thread_ID==1)
			fprintf(stderr, "Process input file: %s\n", filelist1[fi].c_str());
	while ((gzinfile && Read_Tag_gz(Ffilelist1[fi].gzfp,Tail_File.gzfp, Ffilelist1[fi].FILETYPE,R, M)) || (!gzinfile && Read_Tag(Ffilelist1[fi].Input_File,Tail_File.Input_File, Ffilelist1[fi].FILETYPE,R, M)))
	{
		RECOVER_N=0;
		R.Tag_Number=FIRST_READ;
		int SEG_SIZE=75;
		if(SEG_SIZE>=Read_Length) SEG_SIZE=Read_Length-1;
		int SHIFT_SEG_SIZE=(2*SEG_SIZE>Read_Length)? Read_Length-SEG_SIZE-1:SEG_SIZE;
		int Hits_Printed=0;
		Progress++;Total_Reads++;
		Actual_Tag++;
		if (READS_TO_PROCESS && READS_TO_PROCESS <= Total_Reads) break;
		if (ESTIMATE && READS_TO_ESTIMATE <= Total_Reads) break;

		if(gzinfile){
			if ( Progress>=100000 && Thread_ID==1 ) {
				fprintf_time(stderr, "Processed %d reads.\n", Total_Reads);
				Progress = 0;
			}
		}
		else if (Thread_ID==1 && Progress>Number_of_Tags && PROGRESSBAR) 
		{
			off64_t Current_Pos=ftello64(Ffilelist1[fi].Input_File);
			off64_t Average_Length=Current_Pos/Actual_Tag+1;
			Number_of_Tags=(Ffilelist1[fi].File_Size/Average_Length)/100;
			Progress=0;
			Show_Progress(Current_Pos*100/Ffilelist1[fi].File_Size);
		}
		// > 0.15
		if(DO_INDEL_LARGE && Total_Reads>10000 && Bindel>Progress*0.6) DO_INDEL_LARGE=false;

		//Read Head start ..
		R.Real_Len=0;
		for(;R.Tag_Copy[R.Real_Len]!=0 && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
		M.Real_Len=0;
		
		//if(R.Real_Len >= DetectIndelLen) DO_INDEL_LARGE = true;
		//else DO_INDEL_LARGE = false;

                //-------------batmeth1------------------------------------------

                if(R.Real_Len != L_batmeth.STRINGLENGTH)
                {
                        GET_LEN(L_batmeth,R);
                        Split_Read(L_batmeth.STRINGLENGTH_ORG,L_batmeth);//calculate read portions for Batman algo..
                        O_CT.Length_Array[0]=O_CT.Length_Array[1]=O_CT.Length_Array[2]=L_batmeth.STRINGLENGTH;
                        O_RevCT=O_CT;
                        O_GA=O_CT;O_RevGA=O_CT;
                        B_batmeth.IGNOREHEAD=L_batmeth.IGNOREHEAD;B_batmeth.StringLength=L_batmeth.STRINGLENGTH;
                }
                //Read Head start ..

		//-----------------------output seq header --------------------------------moxian
		R.Tag_Copy[R.Real_Len]='\0';
		R.Quality[R.Real_Len]='\0';
		strcpy(R.Raw_Tag_Copy,R.Tag_Copy);
		//
                revR=R_CT=R_GA=R;
                int i;int Lens=R.Real_Len;
                for (i=0;i<=Lens-1;i++)
                {
                        revR.Tag_Copy[Lens-1-i]=Char_To_CharC[R.Tag_Copy[i]];
                };
                Lens = M.Real_Len;
                for (i=0;i<=Lens-1;i++)
                        revM.Tag_Copy[Lens-1-i]=Char_To_CharC[M.Tag_Copy[i]];

                Reverse_Quality(revR.Quality,R,R.Real_Len);
                Reverse_Quality(revM.Quality,M,M.Real_Len);

		FreeQ_Pair(Alignments_Reslut);
		FreeQ_Hits(Align_Hits_CT);FreeQ_Hits(Align_Hits_CT_P);
		Align_Hits_CT.nMisNow = -1;
		//===init
		//printf("\n%s %s %d\n", R.Description, R.Tag_Copy, Thread_ID);
		int align_mis=0;
		char source1 = '1', source2 = '0';
		int cutfortrcalQS = 30;
//------------------------------------------------------------------------------------------------
		//-------------13
		int SecondN = 0;
		int Paired_Score=INT_MAX;
		Mapping_Batmis(align_mis,Paired_Score, source1, source2,R_CT,M_GA,B_batmeth,MF_CT,MC_CT,MF_CT2,MC_CT2,mF_CT,mC_CT,L_batmeth,fwfmiCT,revfmiCT,R,M,Align_Hits_CT,Align_Hits_CT_P,O_CT,O_RevCT,Original_Text_CT,Single_File,Alignments_Reslut,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_CT,RQHALF_CT,RQ_CT);
		source1 = '4', source2 = '0';
		Mapping_Batmis(align_mis,Paired_Score,source1,source2,R_CT,M_GA,B_batmeth,MF_GA,MC_GA,MF_GA2,MC_GA2,mF_GA,mC_GA,L_batmeth,fwfmiGA,revfmiGA,R,M,Align_Hits_CT,Align_Hits_CT_P,O_RevGA,O_GA,Original_Text_GA,Single_File,Alignments_Reslut,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_GA,RQHALF_GA,RQ_GA);
		if(!Alignments_Reslut.empty() && !ESTIMATE )
		{    
			Alignment_Pair pai=Alignments_Reslut.top();
			if((pai.Mismatch==0 && pai.Indel==0) && pai.Clip ==0 && Alignments_Reslut.size()>1 ) 
			{
				ntemp++;   
				Remove_Dup_Final(Alignments_Reslut, SecondN);
				pai=Alignments_Reslut.top();
                                Alignments_Reslut.pop(); Alignment_Pair pai2; pai2.Score = 0, pai2.swscore1 = 0;
                                if(Alignments_Reslut.size()>0) {
					pai2=Alignments_Reslut.top();
					if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
					if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore1 = 0;}
					//if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) 
					{int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
				}
				if(pai.QualityScore>20) passed20++;
				Total_Mapped++;
				if(pai.source1 == '1') 
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				else if(pai.source1 == '4')
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				continue;
			}    
		}
       //////////////////////////////mis 1
       align_mis=1;
		source1 = '1', source2 = '0';
		Mapping_Batmis(align_mis,Paired_Score,source1,source2,R_CT,M_GA,B_batmeth,MF_CT,MC_CT,MF_CT2,MC_CT2,mF_CT,mC_CT,L_batmeth,fwfmiCT,revfmiCT,R,M,Align_Hits_CT,Align_Hits_CT_P,O_CT,O_RevCT,Original_Text_CT,Single_File,Alignments_Reslut,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_CT,RQHALF_CT,RQ_CT);
		source1 = '4', source2 = '0';
		Mapping_Batmis(align_mis,Paired_Score,source1,source2,R_CT,M_GA,B_batmeth,MF_GA,MC_GA,MF_GA2,MC_GA2,mF_GA,mC_GA,L_batmeth,fwfmiGA,revfmiGA,R,M,Align_Hits_CT,Align_Hits_CT_P,O_RevGA,O_GA,Original_Text_GA,Single_File,Alignments_Reslut,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_GA,RQHALF_GA,RQ_GA);
		if(!Alignments_Reslut.empty() && !ESTIMATE )
		{    
			Alignment_Pair pai=Alignments_Reslut.top();
			if(pai.Mismatch<=1 && pai.Indel==0  && pai.Clip ==0 && Alignments_Reslut.size()>1 || (pai.Mismatch==0) )
			{
				ntemp++;   
				Remove_Dup_Final(Alignments_Reslut, SecondN);
				pai=Alignments_Reslut.top();
                                Alignments_Reslut.pop(); Alignment_Pair pai2;pai2.Score = 0, pai2.swscore1 = 0;
                                if(Alignments_Reslut.size()>0) {
					pai2=Alignments_Reslut.top();
					if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
					if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore = 0;}
					//if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) 
					{int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
				}
				if(pai.QualityScore>20) passed20++;
				Total_Mapped++;
				if(pai.source1 == '1')
                    fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
                else if(pai.source1 == '4')
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				continue;
			}
		}
		if(ESTIMATE) continue;
		///////////////////mis 2
		align_mis=2;
		source1 = '1', source2 = '0';
		Mapping_Batmis(align_mis,Paired_Score,source1,source2,R_CT,M_GA,B_batmeth,MF_CT,MC_CT,MF_CT2,MC_CT2,mF_CT,mC_CT,L_batmeth,fwfmiCT,revfmiCT,R,M,Align_Hits_CT,Align_Hits_CT_P,O_CT,O_RevCT,Original_Text_CT,Single_File,Alignments_Reslut,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_CT,RQHALF_CT,RQ_CT);
		source1 = '4', source2 = '0';
		Mapping_Batmis(align_mis,Paired_Score,source1,source2,R_CT,M_GA,B_batmeth,MF_GA,MC_GA,MF_GA2,MC_GA2,mF_GA,mC_GA,L_batmeth,fwfmiGA,revfmiGA,R,M,Align_Hits_CT,Align_Hits_CT_P,O_RevGA,O_GA,Original_Text_GA,Single_File,Alignments_Reslut,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Mishit_File,Pairs,Conversion_Factor,Entries_GA,RQHALF_GA,RQ_GA);
		int mis2N = Align_Hits_CT.AlignTmp.size();
		if(!Alignments_Reslut.empty() && !ESTIMATE )
		{    
			Alignment_Pair pai=Alignments_Reslut.top();
			if(Align_Hits_CT.AlignTmp.size()>1 && pai.Mismatch<=2 && pai.Clip ==0)
			{
				ntemp++;   
				Remove_Dup_Final(Alignments_Reslut, SecondN);
				pai=Alignments_Reslut.top();
                                Alignments_Reslut.pop(); Alignment_Pair pai2; pai2.Score = 0, pai2.swscore1 = 0;
                                if(Alignments_Reslut.size()>0) {
					pai2=Alignments_Reslut.top();
					if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
					if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore1 = 0;}
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) {int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
				}
				if(pai.QualityScore>20) passed20++;
				Total_Mapped++;
				if(pai.source1 == '1')
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				else if(pai.source1 == '4')
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				continue;
			}    
		}
		////////////// indel, mismatch 0
	if(DO_INDEL_SMALL)
	{
                int Paired_Ncutoff = 10000;
                source1 = '1', source2 = '0';
                align_mis=0;
		int indel0N = 0;
                ReplaceCtoT(R_CT);ReplaceGtoA(M_GA);//CTread1 GAread2 CTgenome
                Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_CT,RQ_CT,Original_Text_CT,Entries_CT,fwfmiCT,revfmiCT,R_CT,M_GA,B,Conversion_Factor,mF_CT,mC_CT,MF_CT,MC_CT,MF_CT2,MC_CT2,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);
                //-------------------------------4 2-------------------------------------------------
                source1 = '4', source2 = '0';
		Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_GA,RQ_GA,Original_Text_GA,Entries_GA,fwfmiGA,revfmiGA,R_CT,M_GA,B,Conversion_Factor,mF_GA,mC_GA,MF_GA,MC_GA,MF_GA2,MC_GA2,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);
		Remove_Dup_Top(Align_Hits_CT.AlignTmp, 20);
                if(!Alignments_Reslut.empty() && !ESTIMATE )
                {
//printf("\n%d %d\n", Alignments_Reslut.size(), Align_Hits_CT.AlignTmp.size());
                        Alignment_Pair pai=Alignments_Reslut.top();
			indel0N=Align_Hits_CT.AlignTmp.size();
				if( (indel0N>1 || mis2N>0) && ( (pai.Mismatch<=0 && pai.Indel<=1) || ( pai.Mismatch<=3 && pai.Indel==0 ) ) && pai.Clip ==0 )
				{
					Remove_Dup_Final(Alignments_Reslut, SecondN);
					pai=Alignments_Reslut.top();
					Alignments_Reslut.pop(); Alignment_Pair pai2; pai2.Score = 0, pai2.swscore1 = 0;
					if(Alignments_Reslut.size()>0) {
							pai2=Alignments_Reslut.top();
							if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
							if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
							if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore1 = 0;}
							if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) {int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
					}
					if(pai.QualityScore>20) passed20++;
					Total_Mapped++;
				if(pai.source1 == '1' && pai.Mismatch <= MAX_MISMATCHES)
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				else if(pai.source1 == '4' && pai.Mismatch <= MAX_MISMATCHES)
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
                                continue;
                        }
                }
//	if(DO_INDEL_LARGE){
		////////////// indel, mismatch 1
		source1 = '1', source2 = '0';
		align_mis=1;
		ReplaceCtoT(R_CT);ReplaceGtoA(M_GA);//CTread1 GAread2 CTgenome
		Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_CT,RQ_CT,Original_Text_CT,Entries_CT,fwfmiCT,revfmiCT,R_CT,M_GA,B,Conversion_Factor,mF_CT,mC_CT,MF_CT,MC_CT,MF_CT2,MC_CT2,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);
		//-------------------------------4 2-------------------------------------------------
		source1 = '4', source2 = '0';
		Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_GA,RQ_GA,Original_Text_GA,Entries_GA,fwfmiGA,revfmiGA,R_CT,M_GA,B,Conversion_Factor,mF_GA,mC_GA,MF_GA,MC_GA,MF_GA2,MC_GA2,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);
		Remove_Dup_Top(Align_Hits_CT.AlignTmp, 20);
		if(!Alignments_Reslut.empty() && !ESTIMATE )
		{
			Alignment_Pair pai=Alignments_Reslut.top();
			int indel1N = Align_Hits_CT.AlignTmp.size();
			if( (indel1N > 1 || indel0N>0) && ( (pai.Mismatch<=1 && pai.Indel<=1) || ( pai.Mismatch<=3 && pai.Indel==0 ) ) && pai.Clip ==0)
			{
				Remove_Dup_Final(Alignments_Reslut, SecondN);
				pai=Alignments_Reslut.top();
                Alignments_Reslut.pop(); Alignment_Pair pai2; pai2.Score = 0, pai2.swscore1 = 0;
                if(Alignments_Reslut.size()>0) {
					pai2=Alignments_Reslut.top();
					if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
					if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore1 = 0;}
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) {int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
				}
				if(pai.QualityScore>60) continue;
				if(pai.QualityScore>20) passed20++;
				Total_Mapped++;
				if(pai.source1 == '1' && pai.Mismatch <= MAX_MISMATCHES)
                	fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				else if(pai.source1 == '4' && pai.Mismatch <= MAX_MISMATCHES)
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				continue;
			}
		}
		////////////////Final
	if(DO_INDEL_LARGE){
		source1 = '1', source2 = '0';
		align_mis=2;
		Bindel++;
		ReplaceCtoT(R_CT);ReplaceGtoA(M_GA);//CTread1 GAread2 CTgenome
		Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_CT,RQ_CT,Original_Text_CT,Entries_CT,fwfmiCT,revfmiCT,R_CT,M_GA,B,Conversion_Factor,mF_CT,mC_CT,MF_CT,MC_CT,MF_CT2,MC_CT2,MFLH_CT,MCLH_CT,MFLT_CT,MCLT_CT,MFH_CT,MCH_CT,MFT_CT,MCT_CT, L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);
		//-------------------------------4 2-------------------------------------------------
		source1 = '4', source2 = '0';
		Two_Side_Hit_Finding(align_mis,Paired_Ncutoff,source1,source2,R,M,Read_Length,RQHALF_GA,RQ_GA,Original_Text_GA,Entries_GA,fwfmiGA,revfmiGA,R_CT,M_GA,B,Conversion_Factor,mF_GA,mC_GA,MF_GA,MC_GA,MF_GA2,MC_GA2,MFLH_GA,MCLH_GA,MFLT_GA,MCLT_GA,MFH_GA,MCH_GA,MFT_GA,MCT_GA,L,L_Main,L_Half,L_Third,Actual_Tag,Single_File,Mishit_File,Align_Hits_CT,Align_Hits_CT_P,Pairs,0,SEG_SIZE,SHIFT_SEG_SIZE,Alignments_Reslut);

	}//END FOR LARGE INDEL //not here
	}//END FOR INDEL //not here
		if(Alignments_Reslut.empty() && !ESTIMATE)
		{
			;//fprintf(Single_File,"@\n%s\t%s\t%s\n",R.Description,R.Tag_Copy,R.Quality);
			//fprintf(Single_File,"@\n%s\t%s\t%s\n",M.Description,M.Tag_Copy,M.Quality);
		}
		else if(!ESTIMATE)
		{
//printf("\n-- %d\n", Alignments_Reslut.size());
//			while(Alignments_Reslut.size()>0)//rm
//			{//rm
			Alignment_Pair pai=Alignments_Reslut.top();
			if(pai.swscore1 >= R.Real_Len*0.9 || ((pai.swscore1 >= R.Real_Len*0.8 || pai.swscore1>= R.Real_Len*0.6 && Alignments_Reslut.size()<=2) && R.Real_Len <= 75) || (DO_INDEL_LARGE && ((pai.swscore1 - pai.Clip) > R.Real_Len*0.8 || pai.swscore1 > 100)  ))
			{
				READ PrintR,PrintR2;char source2;
				if( pai.source1=='1' || pai.source1=='4' ) 
				{
						PrintR2=R;
						PrintR=M;
				}else 
				{
						PrintR=R;PrintR2=M;
				}
				Remove_Dup_Final(Alignments_Reslut, SecondN);
				pai=Alignments_Reslut.top();
                                Alignments_Reslut.pop(); Alignment_Pair pai2; pai2.Score = 0; pai2.swscore1 = 0;
                                if(Alignments_Reslut.size()>0) {
					pai2=Alignments_Reslut.top();
                                	if(pai2.swscore >= pai.swscore || pai2.Score <= pai.Score) pai.QualityScore = pai.QualityScore<pai2.QualityScore?pai.QualityScore:pai2.QualityScore;
//printf("\nILLLL %d %d, %d %d\n", pai.swscore1, pai2.swscore1, pai.Score, pai2.Score);
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=5) pai.QualityScore = dist(pai.swscore1, pai2.swscore1);
					if((strcmp(pai.chrom1, pai2.chrom1) == 0 && dist(pai2.Loc1, pai.Loc1)<=100)) {pai2.Score = 0; pai2.swscore1 = 0;}
					if( (strcmp(pai.chrom1, pai2.chrom1) != 0 || dist(pai2.Loc1, pai.Loc1)>100) && dist(pai.swscore1, pai2.swscore1)<=30 && pai.QualityScore>=cutfortrcalQS ) {int QS = Cal_mapQ(SecondN,pai2.swscore1, pai.swscore1, pai.Cigar1.c_str(), pai.Mismatch1, R.Real_Len); if(QS < pai.QualityScore) pai.QualityScore = QS;}
				}
				if(pai.QualityScore>20) passed20++;
				Total_Mapped++;
				if(pai.Loc1>0 && pai.source1 == '1')
                	fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), R.Tag_Copy, R.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
				else if(pai.Loc1>0 && pai.source1 == '4')
					fprintf(Single_File,"%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tSA:i:%d\tSW:i:%d\tSS:i:%d\n", pai.Description1+1, pai.Flag1, pai.chrom1, pai.Loc1, pai.QualityScore, pai.Cigar1.c_str(), revR.Tag_Copy, revR.Quality, pai.Mismatch1, pai.Score, pai2.Score, pai.swscore1, pai2.swscore1);
                                
			}
//			}//rm
//		}//END FOR LARGE INDEL
//		}//END FOR INDEL
		}
	}
	if(fi+1 == filelist1.size()) break;
	pthread_mutex_lock(&flists);
	fi++;
	pthread_mutex_unlock(&flists);
    }
	//--------------------------GA----------------------------
	MF_GA.Left_Mishits.clear();MC_GA.Left_Mishits.clear();
	MF_GA.Right_Mishits.clear();MC_GA.Right_Mishits.clear();
	MF_GA.Mismatches_Backward.clear();MC_GA.Mismatches_Backward.clear();
	MF_GA.Mismatches_Forward.clear();MC_GA.Mismatches_Forward.clear();
	MF_GA.Two_Mismatches_At_End_Forward.clear();MC_GA.Two_Mismatches_At_End_Forward.clear();
	MF_GA.Two_Mismatches_At_End.clear();MC_GA.Two_Mismatches_At_End.clear();
	MF_GA.Possible_20.clear();MC_GA.Possible_20.clear();
	MF_GA.Possible_02.clear();MC_GA.Possible_02.clear();
	//---------------------------CT------------------------------
	MF_CT.Left_Mishits.clear();MC_CT.Left_Mishits.clear();
	MF_CT.Right_Mishits.clear();MC_CT.Right_Mishits.clear();
	MF_CT.Mismatches_Backward.clear();MC_CT.Mismatches_Backward.clear();
	MF_CT.Mismatches_Forward.clear();MC_CT.Mismatches_Forward.clear();
	MF_CT.Two_Mismatches_At_End_Forward.clear();MC_CT.Two_Mismatches_At_End_Forward.clear();
	MF_CT.Two_Mismatches_At_End.clear();MC_CT.Two_Mismatches_At_End.clear();
	MF_CT.Possible_20.clear();MC_CT.Possible_20.clear();
	MF_CT.Possible_02.clear();MC_CT.Possible_02.clear();

}

void Remove_Dup_Final(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut, int & secondN)
{
	int Gap = 100;
	std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> Alignnew_Reslut;
        if(Alignments_Reslut.size()<=1)
        {
                return;
        }
        assert(!Alignments_Reslut.empty());
        int bestscore=0, secondscore = 0;
	secondN = 0;
	Alignment_Pair Top_Aln=Alignments_Reslut.top();
//	Alignnew_Reslut.push(Top_Aln);
        Alignments_Reslut.pop();
        Alignment_Pair Sub_Aln=Alignments_Reslut.top();
	if(( strcmp(Sub_Aln.chrom1, Top_Aln.chrom1) == 0 && dist(Sub_Aln.Loc1,Top_Aln.Loc1) == 0) && (dist(Sub_Aln.Score, Top_Aln.Score) == 0)) {
		Top_Aln.QualityScore = Top_Aln.QualityScore < Sub_Aln.QualityScore ? Top_Aln.QualityScore : Sub_Aln.QualityScore;
		if(secondscore < Top_Aln.swscore1 + 1) secondN ++; 
		if(Top_Aln.swscore1 > bestscore) {
			bestscore = Top_Aln.swscore1;
			secondscore = bestscore;
			secondN = 1;
		}
		Alignnew_Reslut.push(Top_Aln);
	}else {
		Alignnew_Reslut.push(Top_Aln);
		if(secondscore < Top_Aln.swscore1 + 1) secondN ++;
		if(Top_Aln.swscore1 > bestscore) {
                        bestscore = Top_Aln.swscore1;
                        secondscore = bestscore;
			secondN = 1;
                }
	}
	
	while(1) {
		if(Alignments_Reslut.empty()) break;
	        if(!Alignments_Reslut.empty() && ( strcmp(Sub_Aln.chrom1, Top_Aln.chrom1) == 0 && dist(Sub_Aln.Loc1,Top_Aln.Loc1) == 0) && (dist(Sub_Aln.Score, Top_Aln.Score) == 0) )
        	{
	                Alignments_Reslut.pop();
                	Sub_Aln=Alignments_Reslut.top();
        	}else if(Alignments_Reslut.size()>=1){
//			Alignnew_Reslut.push(Sub_Aln);
			Top_Aln = Sub_Aln;
			Alignments_Reslut.pop();
			if(Alignments_Reslut.empty()) break;
			Sub_Aln=Alignments_Reslut.top();
		        if(( strcmp(Sub_Aln.chrom1, Top_Aln.chrom1) == 0 && dist(Sub_Aln.Loc1,Top_Aln.Loc1) == 0) && (dist(Sub_Aln.Score, Top_Aln.Score) == 0)) {
		                Top_Aln.QualityScore = Top_Aln.QualityScore < Sub_Aln.QualityScore ? Top_Aln.QualityScore : Sub_Aln.QualityScore;
                		Alignnew_Reslut.push(Top_Aln);
				if(secondscore < Top_Aln.swscore1 + 1) secondN ++;
				if(Top_Aln.swscore1 > bestscore) {
		                        bestscore = Top_Aln.swscore1;
                		        secondscore = bestscore;
					secondN = 1;
		                }
		        }else {
				Alignnew_Reslut.push(Top_Aln);
				if(secondscore < Top_Aln.swscore1 + 1) secondN ++;
				if(Top_Aln.swscore1 > bestscore) {
                        		bestscore = Top_Aln.swscore1;
		                        secondscore = bestscore;
					secondN = 1;
                		}
			}
		}
			
	}
	Alignments_Reslut = Alignnew_Reslut;
}

void FreeQ_Hits(Align_Hit & Align_Hits)
{
	FreeQ(Align_Hits.H0);
	FreeQ(Align_Hits.H1);
	FreeQ(Align_Hits.H2);
	Align_Hits.AM.clear();
}

void Align_Map_One(READ & R,BATREAD & B,char & source,int & mismatch,LEN & L,BWT* fwfmi,BWT* revfmi,MEMX & MF,MEMX & MC,int & Paired_Score,READ & RawR,Align_Hit & Align_Hits,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments, int moveplus, int moveneg)
{
	OUTPUT O1;
	int Next_Mis1=0,Next_Mis2=0;
	int Last_Mis=-1;
	bool batmeth1=true;//used batmeth1 mode to alignment partial-seq
	std::string readString="";//R.Tag_Copy;// no use
	MF.batmeth1=true;MC.batmeth1=true;

     //   ReplaceCtoT(R);
	Process_Read_bat(R,B,MF,MC);
	bool Whole_Len=false;
	int Hits=0;
	if(source=='4')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L,fwfmi,revfmi,MC);
		Next_Mis1=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source,Genome_Offsets,readString,MC,L,FALSE,revfmi,fwfmi,O1,moveplus, moveneg);
	}else if(source=='1')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L,fwfmi,revfmi,MF);
		Next_Mis1=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source,Genome_Offsets,readString,MF,L,FALSE,revfmi,fwfmi,O1,moveplus, moveneg);
	}
	else if(source=='2')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L,fwfmi,revfmi,MF);
		Next_Mis2=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source,Genome_Offsets,readString,MF,L,FALSE,revfmi,fwfmi,O1,moveplus, moveneg);
	}else if(source=='3')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L,fwfmi,revfmi,MC);
		Next_Mis2=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source,Genome_Offsets,readString,MC,L,FALSE,revfmi,fwfmi,O1,moveplus, moveneg);
	}
}

void Align_Map(READ & R,READ & M,BATREAD & B,char & source1,char & source2,int & mismatch,LEN & L_batmeth,BWT* fwfmi,BWT* revfmi,MEMX & MF,MEMX & MC,int & Paired_Score,READ & RawR,READ & RawM,OUTPUT & O1,OUTPUT & O2,Align_Hit & Align_Hits,Align_Hit & Align_Hits_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments_P)
{
	int Next_Mis1=0,Next_Mis2=0;
	int Last_Mis=-1;bool batmeth1=true;
	std::string readString="";//R.Tag_Copy;// no use
     //   ReplaceCtoT(R);
	if( !Process_Read_bat(R,B,MF,MC) ) return;
	bool Whole_Len=true;
	int Hits=0;

	if(source1=='4')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L_batmeth,fwfmi,revfmi,MC);
		Next_Mis1=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source1,Genome_Offsets,readString,MC,L_batmeth,FALSE,revfmi,fwfmi,O1,0,0);
	}else if(source1=='1')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L_batmeth,fwfmi,revfmi,MF);
		Next_Mis1=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawR,Last_Mis,batmeth1,Align_Hits,Alignments,source1,Genome_Offsets,readString,MF,L_batmeth,FALSE,revfmi,fwfmi,O1,0,0);
	}
	
	Last_Mis=-1;
     //   ReplaceGtoA(M);
        if( !Process_Read_bat(M,B,MF,MC) ) return;
	if(source2=='2')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L_batmeth,fwfmi,revfmi,MF);
		Next_Mis2=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawM,Last_Mis,batmeth1,Align_Hits_P,Alignments_P,source2,Genome_Offsets,readString,MF,L_batmeth,FALSE,revfmi,fwfmi,O2,0,0);
	}else if(source2=='3')
	{
		Hits=Map_Strand(Last_Mis,mismatch,maxhits,L_batmeth,fwfmi,revfmi,MC);
		Next_Mis2=Last_Mis;
		Print_Hits(Whole_Len,Paired_Score,RawM,Last_Mis,batmeth1,Align_Hits_P,Alignments_P,source2,Genome_Offsets,readString,MC,L_batmeth,FALSE,revfmi,fwfmi,O2,0,0);
	}
}

void Detect_result(int mismatch,int Indel,bool & Do_Rescue,bool & realign,char source1,char source2,READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,MEMX & mF,MEMX & mC,BWT* fwfmi,BWT* revfmi,READ & RawRTemp,READ & RawRTemp_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments_P,unsigned char* Original_Text,FILE*  Single_File,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments_P)
{
	bool Max_Pass=true;
	Hit_Info H1,H1_P;
	H1.Loc=H1_P.Loc=UINT_MAX;H1.Status=H1_P.Status=UNMAPPED;
	int Quality_Score1=0,Quality_Score1_P=0;
	int Read_Length=RTemp.Real_Len,Read_Length2=RTemp_P.Real_Len;
	Alignment A1,A1_P;
	int MapQ1=-1,MapQ1_P=-1;
	Good_Alignments=Alignments;
	Get_Basic_MapQ(Good_Alignments,A1,MapQ1);
	if(!(MapQ1 == -1 && MapQ1_P== -1))//if at least one end mapped 
	{
		if(MapQ1!= -1 )//Oneside mapped..
		{
			Mate_Rescue_New(mismatch,Indel,Do_Rescue,realign,Alignments_Reslut,RawRTemp,RawRTemp_P,Original_Text,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,Read_Length,Read_Length2,Alignments,Alignments_P,Good_Alignments,Good_Alignments_P,H1,H1_P,Single_File,Quality_Score1,Quality_Score1_P,A1,MapQ1,mF,mC,fwfmi,revfmi);
		}
		return;
	}
}
void combined3(ALIGNMENT_Q & A,ALIGNMENT_Q & b,ALIGNMENT_Q & c) 
{
	ALIGNMENT_Q B=b,C=c;
	while(B.size()>0)
	{
		Alignment aln=B.top();B.pop();
		A.push(aln);
	}
	while(C.size()>0)
	{
		Alignment aln=C.top();C.pop();
		A.push(aln);
	}
}
void combined2(ALIGNMENT_Q & A,ALIGNMENT_Q & b)
{
	ALIGNMENT_Q B=b;
	while(B.size()>0)
	{
		Alignment aln=B.top();B.pop();
		A.push(aln);
	}
}
void Mapping_Batmis(int & align_mismatch,int & Paired_Score,char source1,char source2,READ & R,READ & M,BATREAD & B,MEMX & MF,MEMX & MC,MEMX & MF2,MEMX & MC2,MEMX & mF,MEMX & mC,LEN & L_batmeth,BWT* fwfmi,BWT* revfmi,READ & RawR1,READ & RawM1,Align_Hit & Align_Hits,Align_Hit & Align_Hits_P,OUTPUT & O1,OUTPUT & O2,unsigned char* Original_Text,FILE*  Single_File,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Main,LEN & L_Half,LEN & L_Third,unsigned & Actual_Tag,FILE* Mishit_File,Pair* & Pairs,unsigned & Conversion_Factor,unsigned Entries,RQINDEX & RQHALF,RQINDEX & RQ)
{
	READ RawR = RawR1, RawM = RawM1;
	MF.batmeth1=true;MC.batmeth1=true;mF.ReadExtend=0;
 	mF.batmeth1=true;mC.batmeth1=true;mC.ReadExtend=0;

	
   Alignment A1,B1,C1;
   Alignment A1_P,C1_P,C2_P;
   int MapQ1=1,MapQ2=1,MapQ3=1;
   int MapQ1_P=1,MapQ2_P=1,MapQ3_P=1;
	A1.Loc=A1_P.Loc=UINT_MAX;

	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_P;

	int Read_Length=R.Real_Len;
	ReplaceCtoT(R);
	
	READ RTemp=R;BATREAD BTemp=B;READ RawRTemp=RawR;

	bool FALSEt=false;
	bool TRUEt=true;
	int Read_Length2=M.Real_Len;
	//Hits=0;
	ReplaceGtoA(M);
	READ RTemp_P=M;BATREAD BTemp_P=B;READ RawRTemp_P=RawM;
	
	READ RTemp2=R;BATREAD BTemp2=B;READ RawRTemp2=RawR;
	READ RTemp2_P=M;BATREAD BTemp2_P=B;READ RawRTemp2_P=RawM;
	int mismatch=0;

        if(align_mismatch>1)
        {
                Alignments=Align_Hits.AlignTmp;
                Alignments_P=Align_Hits_P.AlignTmp;
        }

	if(align_mismatch==0) goto Zero;
	else if(align_mismatch==1) goto One;
	else if(align_mismatch==2) goto Two;
	else exit(0);
	
Zero:
	Align_Map(R,M,B,source1,source2,mismatch,L_batmeth,fwfmi,revfmi,MF,MC,Paired_Score,RawR,RawM,O1,O2,Align_Hits,Align_Hits_P,Alignments,Alignments_P);
        Align_Hits.AlignTmp=Alignments;
	if(source1 == '4' && Align_Hits.H0.size()>1)
	{
		Detect_result(0,0,FALSEt,FALSEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Align_Hits.H0,Align_Hits_P.H0,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0) 
		return;
	}
	return;
	
One:
	mismatch=1;
	Align_Map(R,M,B,source1,source2,mismatch,L_batmeth,fwfmi,revfmi,MF,MC,Paired_Score,RawR,RawM,O1,O2,Align_Hits,Align_Hits_P,Alignments,Alignments_P);
	combined3(Alignments,Align_Hits.H0,Align_Hits.H1);
	Align_Hits.AlignTmp=Alignments;
	if( source1 == '4' && (Align_Hits.H0.size()+Align_Hits.H1.size()>1 || Align_Hits.H0.size()>0 ) ) //&& (Align_Hits_P.H0.size()+Align_Hits_P.H1.size()>1 || Align_Hits_P.H0.size()>0 ))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,FALSEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0) 
		{
			Alignment_Pair C=Alignments_Reslut.top();
			if(C.Mismatch<=1 && C.Indel==0)
				return;
		}
	}

	return;

Two:
	mismatch=2;
	Align_Map(R,M,B,source1,source2,mismatch,L_batmeth,fwfmi,revfmi,MF,MC,Paired_Score,RawR,RawM,O1,O2,Align_Hits,Align_Hits_P,Alignments,Alignments_P);
	combined3(Alignments,Align_Hits.H2, Align_Hits.AlignTmp);
        Align_Hits.AlignTmp=Alignments;
	if(source1 == '4' && (Align_Hits.AM.size()>0 || Align_Hits_P.AM.size()>0))
	{
		if(Align_Hits.AM.size()>0 || Align_Hits_P.AM.size()>0)
		{
			RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
			Detect_result(0,0,FALSEt,FALSEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
			if(Alignments_Reslut.size()>0) 
			{
				Alignment_Pair C=Alignments_Reslut.top();
				if( C.Mismatch<=2 && C.Indel==0)
					return;
			}
		}
		
	}
	return;
}
void Get_Top_Hits(ALIGNMENT_Q & A)
{
	int Ncut=1000;
	ALIGNMENT_Q tmp;
}

void combined5(ALIGNMENT_Q & A,ALIGNMENT_Q & B,ALIGNMENT_Q & C,ALIGNMENT_Q & D,ALIGNMENT_Q & E)
{
	while(B.size()>0)
	{
		Alignment aln=B.top();B.pop();
		A.push(aln);
	}
        while(C.size()>0)
        {                    
		Alignment aln=C.top();C.pop();
                A.push(aln);
        }
        while(D.size()>0)
        {
                Alignment aln=D.top();D.pop();
                A.push(aln);        }
        while(E.size()>0)
        {     
                Alignment aln=E.top();E.pop();
                A.push(aln);
        }

}
void Two_Side_Hit_Finding(int & align_mismatch,int & Paired_Ncutoff,char source1,char source2,READ & RawR1,READ & RawM1,int Read_Length,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,READ & M,BATREAD & B,unsigned & Conversion_Factor,MEMX & mF,MEMX & mC,MEMX & MF,MEMX & MC,MEMX & MF2,MEMX & MC2,MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Main,LEN & L_Half,LEN & L_Third,unsigned & Actual_Tag,FILE*  Single_File,FILE* Mishit_File,Align_Hit & Align_Hits,Align_Hit & Align_Hits_P,Pair* & Pairs,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut)
{
	READ RawR = RawR1, RawM = RawM1;
	MF.batmeth1=false;MC.batmeth1=false;mF.ReadExtend=0;
	mF.batmeth1=false;mC.batmeth1=false;mC.ReadExtend=0;
	Alignment A1;
	Alignment A1_P;
	int MapQ1=1;
	int MapQ1_P=1;
	
	Read_Length=R.Real_Len;
	SW_THRESHOLD=80*Read_Length*match/100;
	SEG_SIZE=75;
	if(SEG_SIZE>=Read_Length) SEG_SIZE=Read_Length-1;
	SHIFT_SEG_SIZE=(2*SEG_SIZE>Read_Length)? Read_Length-SEG_SIZE-1:SEG_SIZE;
	int Read_Length2=M.Real_Len;

	A1.Loc=A1_P.Loc=UINT_MAX;
	Hit_Info H1,H1_P;int Quality_Score1,Quality_Score1_P;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	READ RTemp_P=M;BATREAD BTemp_P=B;
	READ RawRTemp_P=RawM;

	Final_Hit Head_Hit,Mate_Hit;
	L_Main.STRINGLENGTH=SEG_SIZE;
	int Half=SMALLSEED;//50;
	Split_Read(Half,L_Half);
	Split_Read(SEG_SIZE/3,L_Third);
	
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp_P;
	
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_P;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_P;
	H1.Status=H1_P.Status=UNMAPPED;
	READ RTemp2=R;BATREAD BTemp2=B;
	READ RawRTemp2=RawR;
	READ RTemp2_P=M;BATREAD BTemp2_P=B;
	READ RawRTemp2_P=RawM;
	bool PRINT=false;bool FALSEt=false;bool TRUEt=true;
//	if(align_mismatch>0)
	{
		Alignments=Align_Hits.AlignTmp;
		Alignments_P=Align_Hits_P.AlignTmp;
	}
	int mismatch=0;
	if(align_mismatch==0) goto Zero;
	else if(align_mismatch==1) goto One;
	else if(align_mismatch==2) goto Two;
//	else if(align_mismatch==3) goto Quarter;
	else exit(0);
	
Zero:	
	Map_One_SEG_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	if(Alignments_tmp.size()>0) {
		combined2(Alignments,Alignments_tmp);
	}
	if(source1 == '4' && ( Alignments_tmp.size()>0) )
	{
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C=Alignments_Reslut.top();
			if( (C.Indel<=1 && C.Mismatch==0) || (C.Indel==0 && C.Mismatch<=2) && C.Clip ==0 ) {
				Align_Hits.AlignTmp=Alignments;
				return;
			}
		}
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_SEG_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF2,MC2,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	if(Alignments_tmp.size()>0) {
		combined2(Alignments,Alignments_tmp);
	}
	if(source1 == '4' && (Alignments_tmp.size()>0 || Alignments_tmp_P.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C=Alignments_Reslut.top();
			if( (C.Indel<=1 && C.Mismatch==0) || (C.Indel==0 && C.Mismatch<=2) && C.Clip==0 ) {
				Align_Hits.AlignTmp=Alignments;
				return;
			}
		}
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLH,MCLH,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);
	if(source1 == '4' && (Alignments_tmp.size()>0 || Alignments_tmp_P.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C=Alignments_Reslut.top();
			if( (C.Indel<=1 && C.Mismatch==0) || (C.Indel==0 && C.Mismatch<=2) && C.Clip == 0 ) {
				Align_Hits.AlignTmp=Alignments;
				return;
			}
		}
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLT,MCLT,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);
	if( (source1 == '4' || Align_Hits.nMisNow == 0) && (Alignments.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C=Alignments_Reslut.top();
			if((C.Indel<=1 && C.Mismatch==0) || (C.Indel==0 && C.Mismatch<=2) && C.Clip == 0 ) {
				Align_Hits.AlignTmp=Alignments;
				return;
			}
		}
	}
	Align_Hits.AlignTmp=Alignments;
	Align_Hits.nMisNow = 0;
	return;
	
One:
	//////////////////////////////////mis 1
	mismatch=1;
	Map_One_SEG_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	combined2(Alignments,Alignments_tmp);
	if(source1 == '4' && (Alignments_tmp.size()>0 || Alignments_tmp_P.size()>0) )
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C = Alignments_Reslut.top();
			if((C.Indel<=1 && C.Mismatch<=1) ) {
				Align_Hits.AlignTmp=Alignments;
				return;
			}
		}
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_SEG_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF2,MC2,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	combined2(Alignments,Alignments_tmp);
	if(source1 == '4' && (Alignments_tmp.size()>0 || Alignments_tmp_P.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C = Alignments_Reslut.top();
			if((C.Indel<=1 && C.Mismatch<=1) && C.Clip == 0 ) { Align_Hits.AlignTmp=Alignments; return;}
		}
	}

	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLH,MCLH,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);
	if(source1 == '4' && (Alignments_tmp.size()>0 || Alignments_tmp_P.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C = Alignments_Reslut.top();
			if(C.Clip == 0 && ( (C.Indel==0 && C.Mismatch<=3) || (C.Indel<=1 && C.Mismatch<=1)) ) { Align_Hits.AlignTmp=Alignments; return; }
		}
	}

	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLT,MCLT,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);
	if( (source1 == '4' || Align_Hits.nMisNow == 1) && (Alignments.size()>0))
	{
		RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
		Detect_result(0,0,FALSEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
		if(Alignments_Reslut.size()>0)
		{
			Alignment_Pair C = Alignments_Reslut.top();
			if(C.Clip == 0 && ( (C.Indel==0 && C.Mismatch<=3) || (C.Indel<=1 && C.Mismatch<=1)) ) { Align_Hits.AlignTmp=Alignments; return; }
		}
	}

	Align_Hits.AlignTmp=Alignments;
	Align_Hits.nMisNow = 1;
	return;
	
Two:
	////////////////////////////////////mis 2
	mismatch=2;
	Map_One_SEG_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF,MC,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	int Best_Mis=INT_MAX,Best_Indel=INT_MAX;
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if((C.Indel<=1 && C.Mismatch<=2) && C.Clip == 0 ) return;
	}
	combined2(Alignments,Alignments_tmp);
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if( ( (C.Indel==0 && C.Mismatch<=4) || (C.Indel<=1 && C.Mismatch<=2)) && C.Clip == 0 ) return;
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_SEG_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF2,MC2,L,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,0,SEG_SIZE,0);
	combined2(Alignments,Alignments_tmp);
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if( ( (C.Indel==0 && C.Mismatch<=4) || (C.Indel<=1 && C.Mismatch<=2)) && C.Clip == 0 ) return;
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Head(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLH,MCLH,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);//Alignments=Alignments+Align_Hits.SH2;
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if( ( (C.Indel==0 && C.Mismatch<=4) || (C.Indel<=1 && C.Mismatch<=2)) && C.Clip == 0 ) return;
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Half_Tail(Align_Hits,mismatch,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MFLT,MCLT,L_Half,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1);
	combined2(Alignments,Alignments_tmp);//Alignments=Alignments+Align_Hits.ST2;
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if( ( (C.Indel==0 && C.Mismatch<=4) || (C.Indel<=1 && C.Mismatch<=2)) && C.Clip == 0 ) return;
	}
	Align_Hits.AlignTmp=Alignments;
//Quarter:	
	//////////////////////Quarter
	int mis_ID=2; //1;
	RawRTemp = RawR; RawRTemp_P = RawM;
   	Map_One_Quart_Head(mis_ID,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF,MC,MFH,MCH,MFT,MCT,L,L_Third,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,SEG_SIZE,SHIFT_SEG_SIZE,Pairs);
	combined2(Alignments,Alignments_tmp);//Alignments=Alignments+Align_Hits.QH1;
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	if(Alignments_Reslut.size()>0)
	{
		Alignment_Pair C=Alignments_Reslut.top();
		if( ( (C.Indel==0 && C.Mismatch<=4) || (C.Indel<=1 && C.Mismatch<=2)) && C.Clip == 0 ) return;
		Best_Mis=C.Mismatch;
		Best_Indel=C.Indel;
	}
	RawRTemp = RawR; RawRTemp_P = RawM;
	Map_One_Quart_Tail(mis_ID,RawR,source1,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,R,B,Conversion_Factor,MF,MC,MFH,MCH,MFT,MCT,L,L_Third,Actual_Tag,Head_Hit,Mishit_File,Alignments_tmp,Good_Alignments_tmp,PRINT,H1,Quality_Score1,SEG_SIZE,SHIFT_SEG_SIZE,Pairs);
	combined2(Alignments,Alignments_tmp);//Alignments=Alignments+Align_Hits.QH2;
	RTemp=RTemp2;RTemp_P=RTemp2_P;BTemp=BTemp2;BTemp_P=BTemp2_P;
	Detect_result(Best_Mis,Best_Indel,TRUEt,TRUEt,source1,source2,RTemp,RTemp_P,BTemp,BTemp_P,mF,mC,fwfmi,revfmi,RawRTemp,RawRTemp_P,Alignments,Alignments_P,Original_Text,Single_File,Alignments_Reslut,Good_Alignments,Good_Alignments_P);
	Align_Hits.AlignTmp=Alignments;

}

void AlignmentHit(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments, char hitType) {
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> A;
	while(!Alignments.empty()) {
		Alignment aln = Alignments.top();
		aln.hitType = hitType;
		Alignments.pop();
		A.push(aln);
	}
}

int Do_Indel(char source,READ & RawR,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT *fwfmi,BWT *revfmi,MEMX & MFLH,MEMX & MCLH,MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,int StringLength,Pair* & Pairs,READ & R,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,const Hit_Info & H)
{
	Ann_Info A;
	int Err=0,Tot_SW_Scans=Alignments.size();
//	char Temp_Current_TagX[MAXDES];
//	char Temp_Current_TagY[MAXDES];
	char *Forward_Read=MFLH.Current_Tag,*Revcomp_Read=MCLH.Current_Tag;
	char *Forward_Read_raw=MFLH.Current_Tag_raw,*Revcomp_Read_raw=MCLH.Current_Tag_raw;
//	MFLH.Least_Mis=MCLH.Least_Mis=INT_MAX;
//	MFLH.Extend=false;MCLH.Extend=false;
	MFH.Extend=false;MCH.Extend=false;
//	MFLT.Extend=false;MCLT.Extend=false;
	MFT.Extend=false;MCT.Extend=false;
//	MEMX Temp_MFL=MFLT,Temp_MCL=MCLT;
	//------------------------Extend Routines Start ----------------------------

	int Plus_HitsX=0,Minus_HitsX=0;
	int Head_Top_Count,Tail_Top_Count;//# of top hits in H/T
	int Mis_StartX=0,Mis_StartY=0;
	bool Exit_Loop=false;
	int Filter=0;
	int Mis_Penalty=(MODE>=FAST) ? Mis_FA_Score:-mismatch;
	int Match_Bonus=(MODE>=FAST) ? 0:match;
	//------------------------Extend Routines End ----------------------------
	SARange *MFH_Top_Start=MFH.Hit_Array,*MCH_Top_Start=MCH.Hit_Array;
	SARange *MFT_Top_Start=MFT.Hit_Array,*MCT_Top_Start=MCT.Hit_Array;
	SARange *MFH_Subopt_Start,*MCH_Subopt_Start;
	SARange *MFT_Subopt_Start,*MCT_Subopt_Start;
	SARange *MFH_Least_Start,*MCH_Least_Start;
	SARange *MFT_Least_Start,*MCT_Least_Start;
	int Plus_Hits=0,Minus_Hits=0;
	unsigned Conversion_Factor=revfmi->textLength-MFH.L.STRINGLENGTH;
	int Max_MM_GAP_Peng=Max_MM_GAP_Adjust;
	if(HEURISTIC) return Err;
	if(EXACT) Max_MM_GAP_Peng=Max_MM_GAP;

	Nindel++;
	//-----------------Pair Top hits ---------------------------------------
	int Top_MisH=Scan(source,MFH,MCH,Max_MM_GAP_Peng,MCH.L,fwfmi,revfmi,0,Head_Top_Count,INT_MAX);
	if(Top_MisH>=0)
	{
		MFH_Subopt_Start= &MFH.Hit_Array[MFH.Hit_Array_Ptr];MCH_Subopt_Start= &MCH.Hit_Array[MCH.Hit_Array_Ptr];
	}
	int Top_MisT=Scan(source,MFT,MCT,Max_MM_GAP_Peng,MCT.L,fwfmi,revfmi,0,Head_Top_Count,INT_MAX);
	if(Top_MisT>=0)
	{
		MFT_Subopt_Start= &MFT.Hit_Array[MFT.Hit_Array_Ptr];MCT_Subopt_Start= &MCT.Hit_Array[MCH.Hit_Array_Ptr];
	}
NN_indels1+=(MFH.Hit_Array_Ptr+MCH.Hit_Array_Ptr);
NN_indels2+=(MFT.Hit_Array_Ptr+MCT.Hit_Array_Ptr);
	int Current_Penalty=2*(MFT.L.STRINGLENGTH)*Match_Bonus-(Top_MisT+Top_MisH)*Mis_Penalty;
	Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Top_MisT,Top_MisH,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Top_Start,MFT_Top_Start,MCH_Top_Start,MCT_Top_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);
	//-----------------Pair Top hits ---------------------------------------
	if(Alignments.empty() && \
			((Top_MisH==0 && Top_MisT<=Max_MM_GAP_Peng-1)||(Top_MisT==0 && Top_MisH<=Max_MM_GAP_Peng-1)) && \
			((Top_MisT+1<=Max_MM_GAP_Peng) && (Top_MisH+1<=Max_MM_GAP_Peng)) && \
			!SKIPHARD)
	{
		int Sub_MisH=(Top_MisH==-1) ? -1 : Scan(source,MFH,MCH,Max_MM_GAP_Peng,MCH.L,fwfmi,revfmi,Top_MisH+1,Head_Top_Count,INT_MAX);
		if(Sub_MisH>=0)
		{
			MFH_Least_Start= &MFH.Hit_Array[MFH.Hit_Array_Ptr];MCH_Least_Start= &MCH.Hit_Array[MCH.Hit_Array_Ptr];
		}
		int Sub_MisT=(Top_MisT==-1) ? -1 : Scan(source,MFT,MCT,Max_MM_GAP_Peng,MCT.L,fwfmi,revfmi,Top_MisT+1,Head_Top_Count,INT_MAX);
		if(Sub_MisT>=0)
		{
			MFT_Least_Start= &MFT.Hit_Array[MFT.Hit_Array_Ptr];MCT_Least_Start= &MCT.Hit_Array[MCH.Hit_Array_Ptr];
		}
	//-----------------Pair Top-Subopt Start ---------------------------------------
		if(Sub_MisT!= -1)
		{
			Current_Penalty=2*(MFT.L.STRINGLENGTH)*Match_Bonus-(Sub_MisT+Top_MisH)*Mis_Penalty;
			Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Sub_MisT,Top_MisH,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Top_Start,MFT_Subopt_Start,MCH_Top_Start,MCT_Subopt_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);
		}
		if(Sub_MisH!= -1)
		{
			Current_Penalty=2*(MFT.L.STRINGLENGTH)*Match_Bonus-(Top_MisT+Sub_MisH)*Mis_Penalty;
			Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Top_MisT,Sub_MisH,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Subopt_Start,MFT_Top_Start,MCH_Subopt_Start,MCT_Top_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);
		}
	//-----------------Pair Top-Subopt end ---------------------------------------
	//-----------------Pair Subopt-Subopt Start ---------------------------------------
		if((Sub_MisH!= -1 && Sub_MisT!= -1) && (Sub_MisT+Sub_MisH<=Max_MM_GAP_Peng))
		{
			Current_Penalty=2*(MFT.L.STRINGLENGTH)*Match_Bonus-(Sub_MisT+Sub_MisH)*Mis_Penalty;
			Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Sub_MisT,Sub_MisH,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Subopt_Start,MFT_Subopt_Start,MCH_Subopt_Start,MCT_Subopt_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);

		}
	//-----------------Pair Subopt-Subopt end ---------------------------------------
		if((Sub_MisH==Max_MM_GAP_Peng-1 && Top_MisT==0)||(Sub_MisT==Max_MM_GAP_Peng-1 && Top_MisH==0))
		{
			int Least_MisH=(Sub_MisH==-1) ? -1 : Scan(source,MFH,MCH,Max_MM_GAP_Peng,MCH.L,fwfmi,revfmi,Sub_MisH+1,Head_Top_Count,INT_MAX);
			int Least_MisT=(Sub_MisT==-1) ? -1 : Scan(source,MFT,MCT,Max_MM_GAP_Peng,MCT.L,fwfmi,revfmi,Sub_MisT+1,Head_Top_Count,INT_MAX);
			//-----------------Pair Top-Least Start ---------------------------------------
			if(Least_MisH==Max_MM_GAP_Peng && Top_MisT==0)
			{
				Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Top_MisT,Least_MisH,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Least_Start,MFT_Top_Start,MCH_Least_Start,MCT_Top_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);

			}
			if(Least_MisT==Max_MM_GAP_Peng && Top_MisH==0)
			{
				Paired_Extension(Forward_Read_raw,Revcomp_Read_raw,RawR,Original_Text,Entries,revfmi,Top_MisH,Least_MisT,Forward_Read,Revcomp_Read,RQ,Pairs,MFH_Top_Start,MFT_Least_Start,MCH_Top_Start,MCT_Least_Start,StringLength,R,Err,Conversion_Factor,Alignments,Current_Penalty,Tot_SW_Scans);
			}
			//-----------------Pair Top-Least End ---------------------------------------
		}
	}
	AlignmentHit(Alignments, source);
	return Err;
}

#define dist(x,y) (((x)>(y)) ? (x)-(y): (y)-(x))
bool Report_Mismatch_Hits(READ & RawR,unsigned char* Original_Text,READ & R,Final_Hit &  Single_File,const int StringLength,Hit_Info & Mismatch_Hit,int Quality_Score)
{
	assert (Mismatch_Hit.Loc || Mismatch_Hit.Status==UNMAPPED);
	if(Mismatch_Hit.Status==UNIQUEHIT || Mismatch_Hit.Status==SHARP_UNIQUEHIT)//Hit found and uniquely..
	{
		assert(Quality_Score<=40);
		Print_Sam(RawR,Original_Text,Single_File,R,Mismatch_Hit,StringLength,30,Default_Alignment,Mismatch_Hit.Clip_H,Mismatch_Hit.Clip_T,Mismatch_Hit.Cigar);return true;
	}
	else if(Mismatch_Hit.Status==MULTI_HIT)
	{
		Mismatch_Hit.Loc=Mismatch_Hit.Sub_Opt_Hit;
		Print_Sam(RawR,Original_Text,Single_File,R,Mismatch_Hit,StringLength,30,Default_Alignment,Mismatch_Hit.Clip_H,Mismatch_Hit.Clip_T,Mismatch_Hit.Cigar);return true;
	}
	else if(Mismatch_Hit.Status==UNRESOLVED_HIT)
	{
		Mismatch_Hit.Status=MULTI_HIT;
		Print_Sam(RawR,Original_Text,Single_File,R,Mismatch_Hit,StringLength,0,Default_Alignment,0,0,Mismatch_Hit.Cigar);return true;
	}
	else
	{
		assert(Mismatch_Hit.Loc!=UINT_MAX || Mismatch_Hit.Status==UNMAPPED);//||Mismatch_Hit.Status==UNRESOLVED_HIT);
		return false;
	}
	assert(false);
}

bool Report_Single_SW(READ & RawR,char source,unsigned char* Original_Text,const int Err,READ & R,Final_Hit & Printed_Hit,const int StringLength,BATREAD & Read,Hit_Info & Mismatch_Hit,bool & Print_Status,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,int Clip_H,int Clip_T,char *CIG)
{
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	Ann_Info Ann;
	Hit_Info H;
	if(Mismatch_Hit.Status==SW_RECOVERED)
	{
		A=Alignments.top();A.Score= -A.Score;
		if(!CIG)//If cigar not present, multi hit recovery in mismatch state did not actually produce many hits..
		{
			assert(BOOST);//due to clipping of some realignings being stored..
			CIG=A.Cigar;Clip_T=A.Clip_T;Clip_H=A.Clip_H;
		}
		assert(A.Realigned==1);
		H.Org_Loc=A.Loc;H.Sign=A.Sign;H.QScore=A.QScore;H.Status=SW_RECOVERED;
		H.Loc = A.Loc; H.Indel=A.Indel; H.Clip=A.Clip; H.hitType = A.hitType;H.BQScore=A.BQScore;
		Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;

		if (H.Loc+StringLength <= Ann.Size && Err<=1)
		{
			if(H.hitType == '1')
				Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Printed_Hit,R,true,30,A,Clip_H,Clip_T,CIG);
			else if(H.hitType == '4')
				Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Printed_Hit,R,true,30,A,Clip_H,Clip_T,CIG);
			else{
				fprintf(stderr, "\nUne Error 0!\n"); exit(0);
			}
			return true;
		}
		else
		{
			Printed_Hit.CIG="";
			return false;
		}
	}
	else
	{
		int SW_Quality_Score=INT_MAX;
		assert(Alignments.size());
		A=Alignments.top();
		assert(A.Realigned==NO || A.Realigned==1);

		FreeQ(Good_Alignments);
		H.Org_Loc=A.Loc;H.Loc = A.Loc;H.hitType = A.hitType;
		Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;
		H.Sign=A.Sign;H.Mismatch=A.Mismatch;H.Indel=A.Indel;H.Score=A.Score;H.QScore=A.QScore;
		H.Clip=A.Clip;
		if(A.Realigned==NO)/* hit only from indel stage,or 75bp hit..*/
		{
			//assert(Mismatch_Hit.QScore== -1);
			bool Do_Smith_Waterman=true;
			if(H.Sign=='+')
			{
				Hit_Info H2=H;
				Alignment A=RealignFast(RawR,source,Original_Text,H2,StringLength,R,0,0,true);
				if(A.Score!=INT_MAX)
				{
					Do_Smith_Waterman=false;
					assert(A.Score<=0);
					A.Realigned=1;A.Clip_T=A.Clip_H=0;
					A.Extend=1;
					Good_Alignments.push(A);
				}
			}

			if(Do_Smith_Waterman)
			{
				RealignX(RawR,Original_Text,H,Read,StringLength,R,false,Good_Alignments,NULL,Clip_H,Clip_T);
			}
			////RealignX(H,Read,StringLength,R,false,Alignments,Good_Alignments,NULL,Clip_H,Clip_T);

			if(!Good_Alignments.empty())
			{
				A=Good_Alignments.top();
				A.Loc--;
			}
			else
			{
				Printed_Hit.CIG="";
				return false;
			}
		}
		A.Score= -A.Score;

		if(Mismatch_Hit.Status !=UNMAPPED && (A.Score>=Mismatch_Hit.Score))
		{
			Report_Mismatch_Hits(RawR,Original_Text,R,Printed_Hit,StringLength,Mismatch_Hit,30);
			return true;
		}
		else//hit in indel stage..
		{
			assert(A.QScore!=INT_MAX);
			H.Org_Loc=A.Loc;H.Loc = A.Loc;H.Sign=A.Sign;H.QScore=A.QScore;H.Status=UNIQUEHIT;
			H.hitType = A.hitType;H.BQScore=A.BQScore;
			//A.Score= -A.Score;
			Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;
			assert(A.Realigned);
			if (H.Loc+StringLength <= Ann.Size && Err<=1)
			{
				if(H.hitType == '1')
					Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Printed_Hit,R,true,30,A,A.Clip_H,A.Clip_T,A.Cigar);
				else if(H.hitType == '4')
					Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Printed_Hit,R,true,30,A,A.Clip_H,A.Clip_T,A.Cigar);
				else {
					fprintf(stderr, "\nUne Error 1 !\n");
					exit(0);
				}
				return true;
			}
			else
			{
				Printed_Hit.CIG="";
				return false;
			}
		}
	}
	assert(false);
	//return true;
}

void Get_Best_Alignment_Pair(READ & RawR,char source,unsigned char* Original_Text,Alignment & A,Alignment & B,READ & R,const int StringLength,BATREAD & Read,Hit_Info & H,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,int & Clip_H,int & Clip_T,char* CIG,bool PRINT,bool DUMMY_FORCED)
{
	int C=0;

	FreeQ(Good_Alignments);
	int Best=INT_MAX;
	int Inc_Count=0;
	char CIG2[MAX_SIGLEN];
	int TClip_T,TClip_H;
	//int Clip_T,Clip_H;
	while(!Alignments.empty())
	{

		A=Alignments.top();Alignments.pop();
		assert(A.Realigned==NO || A.Realigned==1);
		//if(Force_Indel && !A.Indel && A.QScore!= -1) continue;//no indels or clippings..
		if(MODE!=VERYSENSITIVE && Force_Indel && -A.Score>Best && A.QScore != -1) continue;
		if(BOOST && A.Loc==INT_MAX) continue;

		C++;
		if(C>CUT_MAX_SWALIGN) break;//100000) break;//60) break;
		H.Org_Loc=A.Loc;H.Sign=A.Sign;
		H.hitType = A.hitType;
		if(A.Realigned==NO)//only realign if necessary..
		{
			bool Do_Smith_Waterman=true;
			if(H.Sign=='+')
			{
				Hit_Info H2=H;
				A=RealignFast(RawR,source,Original_Text,H2,R.Real_Len,R,0,0,true);
				if(A.Score!=INT_MAX)
				{
					Do_Smith_Waterman=false;
					assert(A.Score<=0);
					A.Realigned=1;A.Clip_T=A.Clip_H=0;
					A.Extend=1;
					//Good_Alignments.push(A);
				}
			}
			if(Do_Smith_Waterman)
			{
				A=RealignX(RawR,Original_Text,H,Read,StringLength,R,true,Good_Alignments,CIG2,TClip_H,TClip_T);//dont push to Q..
			}
		}
		else
		{
			assert(A.Cigar);
			CIG2[0]=0;
		}

		if(A.Sign == '-')
		{
			A.Loc+=(R.Real_Len-StringLength);
		}
		if(A.Loc>1) Good_Alignments.push(A);

		if(A.Loc!=INT_MAX && -A.Score<Best) 
		{
			C=0;
			Best= -A.Score;
			//strcpy(CIG,CIG2);Clip_T=TClip_T;Clip_H=TClip_H;
			strcpy(CIG,A.Cigar);Clip_T=A.Clip_T;Clip_H=A.Clip_H;
		}

	}

	if(!PRINT) return;

	if(!Good_Alignments.empty())//Pick the top two..
	{
		A=Good_Alignments.top();Good_Alignments.pop();A.Score= -A.Score;
		if(!Good_Alignments.empty())
		{
			B=Good_Alignments.top();Good_Alignments.pop();B.Score= -B.Score;
			int Ahl = 0, Atl = 0, Bhl = 0, Btl = 0;
			get_len_match(A, Ahl, Atl);
			get_len_match(B, Bhl, Btl);
			while( (Align_comp(A, B, Ahl, Atl, Bhl, Btl, R.Real_Len) || (dist(A.Loc,B.Loc)<10) && !DUMMY_FORCED) )//std::max(10,(R.Real_Len-StringLength))) 
			{
				if(Good_Alignments.empty())
				{
					B.Score=INT_MAX; 
					break;
				}
				else
				{
					B=Good_Alignments.top();Good_Alignments.pop();B.Score= -B.Score;
					get_len_match(B, Bhl, Btl);
				}
			}
			if(dist(A.Loc,B.Loc)<10 )//&& !DUMMY_FORCED)//std::max(10,(R.Real_Len-StringLength))) 
				B.Score=INT_MAX; 
		}
		else 
			B.Score=INT_MAX; 
	}
	else A.Score=B.Score=INT_MAX;
	if(CIG[0]!=0)//will not pass thru second alignment..
	{
	//	A.Loc--; //changed by moxian, because is no need to to this
		if(A.Sign=='-')
			A.Loc-=(R.Real_Len-StringLength);
	}
}

int Calc_SW_Quality(READ & RawR,unsigned char* Original_Text,Alignment A,READ & R,const int StringLength,BATREAD & Read,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments)
{
	return 30;
	Alignment B;
	B.Clip_H=B.Clip_T=INT_MAX;
	assert(A.QScore!=INT_MAX);
	//float Top_Prob=pow(10,-A.QScore/10);
	float Top_Prob=Pow10(A.QScore);
	float Prob_Sum=0;
	Hit_Info H;

	if(A.Score>B.Score-10)//use only best hits..
	{
		//Prob_Sum=Top_Prob+pow(10,-B.QScore/10);
		Prob_Sum=Top_Prob+Pow10(B.QScore);
		while(!Good_Alignments.empty())
		{
			B=Good_Alignments.top();Good_Alignments.pop();B.Score= -B.Score;
			if(A.Score<B.Score-10)
			{
				if(B.QScore==INT_MAX)
				{
					H.Org_Loc=B.Loc;H.Sign=B.Sign;H.hitType = B.hitType;
					B=Realign(RawR,Original_Text,H,Read,StringLength,R,true,Good_Alignments);
					B.Score= -B.Score;
				}
				assert(B.QScore!=INT_MAX);
				//Prob_Sum+=pow(10,-B.QScore/10);
				Prob_Sum+=Pow10(B.QScore);
			}
			else
				break;
		}
		float Mapping_Quality,Posterior=Top_Prob/(Prob_Sum);
		Mapping_Quality=(int(Posterior)==1) ? 30 : -10*log10(1-Posterior);
		assert(int(Mapping_Quality)>=0);
		return std::max(1,int(Mapping_Quality));
	}
	else
	{
		return 30;
	}

}

bool Report_SW_Hits(READ & RawR,char source,unsigned char* Original_Text,const int Err,READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,Hit_Info & Mismatch_Hit,int Quality_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool Force_Indel,bool PRINT,bool DUMMY_FORCED)
{
	int Alignment_Count=Alignments.size();
	bool Print_Status=false;
	if(!Alignment_Count)//Smith waterman recovery failed or not performed..(Mismatch stage)
	{
		if(!PRINT) 
		{
			return true;
		}
		Print_Status=Report_Mismatch_Hits(RawR,Original_Text,R,Single_File,StringLength,Mismatch_Hit,Quality_Score);
		//Print_Status=true;
	}
	else
	{
		Alignment A,B;
		Ann_Info Ann;
		Hit_Info H;H.Sub_Opt_Score=INT_MAX;H.BQScore=INT_MAX;

		int SW_Quality_Score=INT_MAX;
		int Flag;
		int Clip_H,Clip_T;
		char CIG[MAX_SIGLEN];

		if(Alignment_Count==1)//Smith waterman got a unique hit..
		{
			if(DUMMY_FORCED)
			{
				Alignment A=Alignments.top();A.Score=-A.Score;
				if(A.Loc==INT_MAX) return false;
				A.Sub_Opt_Score=H.Sub_Opt_Score=A.Score+21;
				strcpy(CIG,A.Cigar);Clip_T=A.Clip_T;Clip_H=A.Clip_H;
				
				int SW_Quality_Score=Calc_SW_Quality(RawR,Original_Text,A,R,StringLength,Read,Alignments,Good_Alignments);
				H.Org_Loc=A.Loc;H.Loc = A.Loc;Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;H.Sign=A.Sign;H.Mismatch=A.Mismatch;H.Indel=A.Indel;
				H.Clip=A.Clip;
				H.Score=A.Score;H.QScore=A.QScore;H.SW_Sub_Opt_Score=A.SW_Score-20;
				H.BQScore=A.BQScore;
				H.swscore = A.swscore; H.hitType = A.hitType;
				if(Err>1) Flag=4; else if (H.Sign=='+') Flag=0; else Flag=16; 
				if(TOP_TEN) if (H.Sign=='+') Flag=0; else Flag=16; 
				if(Flag!=4||TOP_TEN) 
				{
					if(Mismatch_Hit.Status==PAIRED_SW)
						H.Status=Mismatch_Hit.Status;
					else
						H.Status=Mismatch_Hit.Status=MULTI_HIT;
					if(H.hitType == '1')
						Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Single_File,R,true,SW_Quality_Score,A,Clip_H,Clip_T,CIG);
					else if(H.hitType == '4')
						Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Single_File,R,true,SW_Quality_Score,A,Clip_H,Clip_T,CIG);
					else {
						fprintf(stderr, "\nUne Error 2 !\n"); exit(0);
					}
					Mismatch_Hit = H;
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				if(!PRINT)
				{
					Get_Best_Alignment_Pair(RawR,source,Original_Text,A,B,R,StringLength,Read,H,Alignments,Good_Alignments,Force_Indel,Clip_H,Clip_T,CIG,PRINT,DUMMY_FORCED);
					Mismatch_Hit = H;
					return true;
				}
				if(!Report_Single_SW(RawR,source,Original_Text,Err,R,Single_File,StringLength,Read,Mismatch_Hit,Print_Status,Alignments,Good_Alignments,0,0,NULL)) 
					return false;
				else 
					return true;
			}
		}
		else if (Alignment_Count)//Multiple SW hits..
		{	
			Get_Best_Alignment_Pair(RawR,source,Original_Text,A,B,R,StringLength,Read,H,Alignments,Good_Alignments,Force_Indel,Clip_H,Clip_T,CIG,PRINT,DUMMY_FORCED);
//		A=Alignments.top();Alignments.pop();B=Alignments.top();FreeQ(Alignments);
//	if(A.Score<0) A.Score=-A.Score;
//	if(B.Score<0) B.Score=-B.Score;
			if(!PRINT) return true;
//printf("\neee %ld %d %ld %d %d %d %s %s %d %d\n",A.Loc,A.Score,B.Loc,B.Score,A.QScore,B.QScore, A.Cigar, B.Cigar, A.swscore, B.swscore);
//			if(A.swscore > B.swscore + 30) B.Score=INT_MAX;
//			if(A.Mismatch==B.Mismatch)
			{
				if(dist(A.Score,B.Score)<10)//Same mismatch, too close to call..
				{
					A.Score=B.Score;
				}
			}
			if(Mismatch_Hit.Status!=SW_RECOVERED && Mismatch_Hit.Status!=UNMAPPED && A.Score>=Mismatch_Hit.Score && Mismatch_Hit.Status!=PAIRED_SW)
			{
				assert(false);
				Print_Status=Report_Mismatch_Hits(RawR,Original_Text,R,Single_File,StringLength,Mismatch_Hit,30);
			}
			else if(B.Score==INT_MAX)
			{
				A.Score= -A.Score;A.Sub_Opt_Score=INT_MAX; Alignments.push(A);
				assert(A.Score!=INT_MAX);
				if(!Report_Single_SW(RawR,source,Original_Text,Err,R,Single_File,StringLength,Read,Mismatch_Hit,Print_Status,Alignments,Good_Alignments,Clip_H,Clip_T,CIG)) 
					return false;
				else
					return true;
			}
			else if(A.Score<B.Score)//Top hit with multiple hits..
			{
				if(A.QScore==INT_MAX)
				{
					H.Org_Loc=A.Loc;H.Sign=A.Sign;H.hitType = A.hitType;
					A=Realign(RawR,Original_Text,H,Read,StringLength,R,true,Good_Alignments);
					A.Score= -A.Score;
				}
				A.Sub_Opt_Score=B.Score;
				SW_Quality_Score=Calc_SW_Quality(RawR,Original_Text,A,R,StringLength,Read,Alignments,Good_Alignments);
				H.Org_Loc=A.Loc;H.Loc = A.Loc;Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;H.Sign=A.Sign;H.Mismatch=A.Mismatch;H.Indel=A.Indel;
				H.Clip=A.Clip;
				H.Score=A.Score;H.Sub_Opt_Score= B.Score;H.QScore=A.QScore;H.SW_Sub_Opt_Score=B.SW_Score;
				H.swscore = A.swscore; H.hitType = A.hitType;
				if(Err>1) Flag=4; else if (H.Sign=='+') Flag=0; else Flag=16; 
				if(TOP_TEN) if (H.Sign=='+') Flag=0; else Flag=16; 
			//strcpy(CIG,A.Cigar);Clip_T=A.Clip_T;Clip_H=A.Clip_H;
				if(Flag!=4||TOP_TEN) 
				{
					H.Status=Mismatch_Hit.Status=MULTI_HIT;
					if(H.hitType == '1')
						Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Single_File,R,true,SW_Quality_Score,A,Clip_H,Clip_T,CIG);
					else if(H.hitType == '4')
						Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Single_File,R,true,SW_Quality_Score,A,Clip_H,Clip_T,CIG);
					else {
						fprintf(stderr, "\nUne Error 3 !\n");
						exit(0);
					}
					Mismatch_Hit = H;
					Print_Status=true;
				}
				if(Flag==4)
				{
					return false;
				}
				if(TOP_TEN)
				{
					Report_Single(RawR,Original_Text,R,Single_File,StringLength,Read,Print_Status,Clip_H,Clip_T,B);
					int c=0;std::set <unsigned> S;
					S.insert(A.Loc);S.insert(B.Loc);
					while(TOP_TEN >2 && !Good_Alignments.empty() && c!=(TOP_TEN-2))
					{
						A=Good_Alignments.top();Good_Alignments.pop();A.Score= -A.Score;
						if(S.find(A.Loc)!=S.end())
							continue;
						S.insert(A.Loc);
						Report_Single(RawR,Original_Text,R,Single_File,StringLength,Read,Print_Status,Clip_H,Clip_T,A);
						c++;
					}
				}
			}
			else //Multiple hits..
			{
				H.Org_Loc=A.Loc;H.Loc = A.Loc;Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;H.Sign=A.Sign;H.Mismatch=A.Mismatch;H.Indel=A.Indel;H.Score=A.Score;H.QScore=A.QScore;H.swscore=A.swscore;H.BQScore=A.BQScore;
				H.Clip=A.Clip; H.hitType = A.hitType;
				if (H.Sign=='+') Flag=0; else Flag=16; 
				if(H.hitType == '1')
					Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Single_File,R,true,0,A,Clip_H,Clip_T,CIG);
				else if(H.hitType == '4')
					Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Single_File,R,true,0,A,Clip_H,Clip_T,CIG);
				else {
					fprintf(stderr, "\nUne Error 4!\n"); exit(0);
				}
				Mismatch_Hit = H;
				Print_Status=true;
				if(TOP_TEN)
				{
					Report_Single(RawR,Original_Text,R,Single_File,StringLength,Read,Print_Status,Clip_H,Clip_T,B);
					int c=0;std::set <unsigned> S;
					S.insert(A.Loc);S.insert(B.Loc);
					while(TOP_TEN >2 && !Good_Alignments.empty() && c!=(TOP_TEN-2))
					{
						A=Good_Alignments.top();Good_Alignments.pop();A.Score= -A.Score;
						if(S.find(A.Loc)!=S.end())
							continue;
						S.insert(A.Loc);
						Report_Single(RawR,Original_Text,R,Single_File,StringLength,Read,Print_Status,Clip_H,Clip_T,A);
						c++;
					}
				}
			}

		}
		else assert(false); 
	}
	return Print_Status;
}

void  Paired_Extension(char* Fwd_Read_raw,char *Revcomp_Read_raw,READ & Raw,unsigned char* Original_Text,unsigned Entries,BWT *revfmi,int Last_MisT,int Last_MisH,char* Fwd_Read,char *Revcomp_Read, RQINDEX & RQ,Pair* & Pairs,SARange* & MFH_Hit_Array,SARange* & MFT_Hit_Array,SARange* & MCH_Hit_Array,SARange* & MCT_Hit_Array,int StringLength,READ & R,int & Err,unsigned Conversion_Factor,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Penalty,int & Tot_SW_Scans)
{
	if(Last_MisT>=0 && Last_MisH>=0)
	{
		char In_Large;
		Ann_Info A;
		int Pairs_Index=0,Pairings=0;
		Head_Tail(revfmi,MFH_Hit_Array,MFT_Hit_Array,2*StringLength/2+INDELGAP,MAXCOUNT,In_Large,RQ,Entries,Pairs,Pairs_Index,Pairings,Err,Conversion_Factor);
		int SW_Hits=0;
		SW_Hits+=Do_SW_Pair(Fwd_Read_raw+StringLength/3,Raw,Original_Text,Pairs,Fwd_Read+StringLength/3,StringLength/3+INDELGAP,StringLength, Err,StringLength/3+1,'+',Alignments,Current_Penalty,Tot_SW_Scans);
		Pairs_Index=0,Pairings=0;
		Head_Tail(revfmi,MCT_Hit_Array,MCH_Hit_Array,2*StringLength/2+INDELGAP,MAXCOUNT,In_Large,RQ,Entries,Pairs,Pairs_Index,Pairings,Err,Conversion_Factor);
		SW_Hits=0;
		SW_Hits+=Do_SW_Pair(Revcomp_Read_raw+StringLength/3,Raw,Original_Text,Pairs,Revcomp_Read+StringLength/3,StringLength/3+INDELGAP,StringLength, Err,StringLength/3+1,'-',Alignments,Current_Penalty,Tot_SW_Scans);
	}
}

void Extend_Right(char* Temp_Current_Tag_raw,READ & Raw,RQINDEX & RQHALF,unsigned char* Original_Text,int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter )
{
	if(Plus_Hits+Minus_Hits)
	{
		bool Head_Disc=false,Tail_Disc=false;

		if(Plus_Hits)
		{
			int Aln_Index=0;
			int SW_Hits=0;
			SW_Hits+=Do_SW(MFL.Current_Tag_raw+MFL.L.STRINGLENGTH,Raw,RQHALF,Original_Text,revfmi,MFL.Hit_Array,MFL.Current_Tag+MFL.L.STRINGLENGTH,MFL.L.STRINGLENGTH+INDELGAP,StringLength, Err,MFL.L.STRINGLENGTH+1,'+',Current_Score,Alignments,Tot_SW_Scans,Filter);

		}
		if(Minus_Hits)
		{
			int Aln_Index=0;
			int SW_Hits=0;
			SW_Hits+=Do_SW(Temp_Current_Tag_raw+MCL.L.STRINGLENGTH,Raw,RQHALF,Original_Text,revfmi,MCL.Hit_Array,Temp_Current_Tag+MCL.L.STRINGLENGTH,MCL.L.STRINGLENGTH+INDELGAP,StringLength, Err,-MCL.L.STRINGLENGTH-INDELGAP,'-',Current_Score,Alignments,Tot_SW_Scans,Filter);
		}

	}
}

void Extend_Left(char* Temp_Current_Tag_raw,READ & RawR,RQINDEX & RQHALF,unsigned char* Original_Text,int Plus_Hits,int Minus_Hits,BWT* revfmi,MEMX & MFL,MEMX & MCL,char* Temp_Current_Tag,int StringLength,int & Err, READ & R,int Mis_In_Anchor,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Tot_SW_Scans,int & Filter)
{
	if(Plus_Hits+Minus_Hits)
	{
		if(Plus_Hits)
		{
			int Aln_Index=0;
			int SW_Hits=0;
			SW_Hits+=Do_SW(Temp_Current_Tag_raw+MFL.L.STRINGLENGTH,RawR,RQHALF,Original_Text,revfmi,MFL.Hit_Array,Temp_Current_Tag+MFL.L.STRINGLENGTH,MFL.L.STRINGLENGTH+INDELGAP,StringLength,Err,-MFL.L.STRINGLENGTH-INDELGAP,'+',Current_Score,Alignments,Tot_SW_Scans,Filter);
		}
		if(Minus_Hits)
		{
			int Aln_Index=0;
			int SW_Hits=0;
			SW_Hits+=Do_SW(MCL.Current_Tag_raw+MCL.L.STRINGLENGTH,RawR,RQHALF,Original_Text,revfmi,MCL.Hit_Array,MCL.Current_Tag+MCL.L.STRINGLENGTH,MCL.L.STRINGLENGTH+INDELGAP,StringLength, Err,MCL.L.STRINGLENGTH+1,'-',Current_Score,Alignments,Tot_SW_Scans,Filter);
		}
	}
}


//map reads..
int Head_Tail(BWT* revfmi,SARange* Head_Hits,SARange* Tail_Hits,int Insert,int MAXCOUNT,char & In_Large,RQINDEX & R,unsigned Entries,Pair* Pairs,int & Pairs_Index,int & HITS,int & Err,unsigned Conversion_Factor)
{
	assert (MAXCOUNT >0);assert (revfmi);assert(Head_Hits);assert(Tail_Hits);assert(Entries>0);assert(Pairs);
	HITS=0;
	SARange Head,Tail;
	int TGap,HGap;
	if (Pairs_Index >=MAXCOUNT){if(!Err) Err=1;return HITS;}
	for(int i=0;Head_Hits[i].Start && HITS <MAXCOUNT;i++)//Iterate Head Hits
	{
		for(int j=0;Tail_Hits[j].Start && HITS <MAXCOUNT;j++)//With Tail Hits
		{
			Head.Start=Head_Hits[i].Start;
			Head.End=Head_Hits[i].End;
			Head.Mismatches=Head_Hits[i].Mismatches;
			Tail.Start=Tail_Hits[j].Start;
			Tail.End=Tail_Hits[j].End;
			Tail.Mismatches=Tail_Hits[j].Mismatches;
			TGap=Tail.End-Tail.Start;
			HGap=Head.End-Head.Start;
			if(HGap>=INDEX_RESOLUTION || TGap>=INDEX_RESOLUTION){if(!Err) Err=1;}

			if (TGap> MAXGAP || HGap > MAXGAP)//sa ranges not cached..
			{
				In_Large=TRUE;
				continue;
			}
			if (TGap<= SAGAP_CUTOFF || HGap <= SAGAP_CUTOFF)//Small sa ranges
			{
				if(TGap && HGap)
				{
					unsigned Start;
					if(HGap<TGap)
					{
						Start=Head.Start;
						for (int i=0;i<=HGap;i++)
						{
							Head.Start=Start;
							Head.End=Head.Start=Conversion_Factor-BWTSaValue(revfmi,Head.Start);
							Search_Small_Gap(revfmi,R,Head,Tail,Insert,Pairs,Pairs_Index,MAXCOUNT,HITS,Entries,Conversion_Factor);
							if (HITS >= MAXCOUNT) break;
							Start++;
						}
					}
					else
					{
						Start=Tail.Start;
						for (int i=0;i<=TGap;i++)
						{
							Tail.Start=Start;
							Tail.End=Tail.Start=Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
							Search_Small_Gap(revfmi,R,Head,Tail,Insert,Pairs,Pairs_Index,MAXCOUNT,HITS,Entries,Conversion_Factor);
							if (HITS >= MAXCOUNT) break;
							Start++;
						}
					}
				}
				else
				{
					Search_Small_Gap(revfmi,R,Head,Tail,Insert,Pairs,Pairs_Index,MAXCOUNT,HITS,Entries,Conversion_Factor);
				}
			}
			else//Two big SA gaps...
			{
				Get_Head_Tail(revfmi,R,Head,Tail,Insert,Pairs,Pairs_Index,MAXCOUNT,HITS,Entries,Conversion_Factor);
			}
			if (Pairs_Index >=MAXCOUNT){if(!Err) Err=1;return HITS;}

		}
	}
	Pairs[Pairs_Index].Head=0;
	assert (HITS >=0);
	return HITS;
}



void Verbose(char *BWTFILE,char *OCCFILE,char *REVBWTINDEX,char *REVOCCFILE,char *PATTERNFILE,char *HITSFILE,char* LOCATIONFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID)
{
	static int Pass=0;

	if (!MISC_VERB) return;
	if(!Pass)
	{
		fprintf(stderr,"BAT ALIGN - Long ..\n");
		fprintf(stderr,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
		fprintf(stderr,"Using the genome files\n %s\t %s\n %s\t %s\n", BWTFILE,OCCFILE,REVBWTINDEX,REVOCCFILE); 
		fprintf(stderr,"Location file: %s\n", LOCATIONFILE); 
/*		fprintf(stderr,"Query File : %s \t\t Output file: %s\n",PATTERNFILE,HITSFILE);
		if(Patternfile_Count) fprintf(stderr,"Mate File : %s \n ",PATTERNFILE1);
		fprintf(stderr,"Length of Tags: %d\t", L.STRINGLENGTH);
		fprintf(stderr,"Mismatches allowed : %d\n",MAX_MISMATCHES);
		if (FILETYPE == FQ) fprintf(stderr,"FASTQ file..\n"); else fprintf(stderr,"FASTA file..\n");
		if (UNIQ_HIT) fprintf(stderr,"Unique hit mode..\n");
		if (PRIORITYMODE) fprintf(stderr,"Lowest mismatch mode..\n");
		fprintf(stderr,"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n");
*/	}
	Pass++;
}
void FileInfo(char *PATTERNFILE,char *HITSFILE,char MAX_MISMATCHES,int Patternfile_Count,char* PATTERNFILE1,char FILETYPE,LEN & L,char FORCESOLID)
{
	static int Pass2=0;
	if(!Pass2)
	{
		fprintf(stderr,"==============================================================\n");
		fprintf(stderr,"Query File : %s \t\t Output file: %s\n",PATTERNFILE,HITSFILE);
		if(Patternfile_Count) fprintf(stderr,"Mate  File : %s\n ",PATTERNFILE1);
		fprintf(stderr,"Length of Tags: %d\n", L.STRINGLENGTH);
		//fprintf(stderr,"Mismatches allowed : %d\n",MAX_MISMATCHES);
		if (SOLID) 
		{
			if (FORCESOLID) fprintf (stderr,"DIBASE-SOLiD reads...\n");
			else fprintf (stderr,"SOLiD reads...\n");
		}
		if (FILETYPE == FQ) fprintf(stderr,"FASTQ file..\n"); else fprintf(stderr,"FASTA file..\n");
		if (UNIQ_HIT) fprintf(stderr,"Unique hit mode..\n");
		if (PRIORITYMODE) fprintf(stderr,"Lowest mismatch mode..\n");
		fprintf(stderr,"----------------------------------------------------------------------\n");
	}
	Pass2++;
}
void Init(In_File & IN,FMFILES F,RQINDEX R,BATPARAMETERS & BP,char Solid,char Npolicy,LEN & L)
{
	JUMP=BP.INDELSIZE;
	INSERT=DELETE=BP.INDELSIZE;
	INDELGAP=2*JUMP+1;

	SOLID=Solid; 
	srand(0);
	if (SOLID) 
	{
		L.IGNOREHEAD=2;
		PAIRING_SCHEME = 5353;
	}
	//NPOLICY=Npolicy;
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;
	Char_To_Code['0']=0;Char_To_Code['1']=1;Char_To_Code['2']=2;Char_To_Code['3']=3;Char_To_Code['4']=0;
	Char_To_CodeC['A']=3;Char_To_CodeC['C']=2;Char_To_CodeC['G']=1;Char_To_CodeC['T']=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;
		Char_To_CodeC['n']=0;Char_To_CodeC['N']=0;//0 moxian//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC['.']=3;//3; moxian
		Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;
			Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['0']=3;Char_To_CodeC['1']=2;Char_To_CodeC['2']=1;Char_To_CodeC['3']=0;

	Char_To_CharC['A']='T';Char_To_CharC['C']='G';Char_To_CharC['G']='C';Char_To_CharC['T']='A';Char_To_CharC['a']='t';Char_To_CharC['c']='g';Char_To_CharC['g']='c';Char_To_CharC['t']='a';
		Char_To_CharC['n']='n';Char_To_CharC['N']='N';//we are using character count to store the fmicode for acgt
	//Char_To_CodeCS['3']=3;Char_To_CodeCS['2']=2;Char_To_CodeCS['1']=1;Char_To_CodeCS['0']=0;Char_To_CodeCS['.']=0;
	Code_To_CodeC[0]=3;Code_To_CodeC[1]=2;Code_To_CodeC[2]=1;Code_To_CodeC[3]=0;
	//if (SOLID) sprintf(Default_Cigar,"%dM",IN.STRINGLENGTH-1);else sprintf(Default_Cigar,"%dM",IN.STRINGLENGTH);
	if (BP.STD)//set a margin of avg+3*std
	{
		BP.FLANKSIZE= 3*BP.STD;//for normal pairing
		BP.SW_FLANKSIZE= 2*IN.STRINGLENGTH+6*BP.STD;
	}
	else if (!BP.SW_FLANKSIZE)
	{
		BP.SW_FLANKSIZE= IN.STRINGLENGTH/2+50;//check 50 bp either side if not specified...
		if (BP.SW_FLANKSIZE > BP.INSERTSIZE) BP.SW_FLANKSIZE=BP.INSERTSIZE;
	}
	if (BP.SW_FLANKSIZE>=SW_STRING_BUFFER) {fprintf(stderr,"SW Flank too large.. Defaulting\n");BP.SW_FLANKSIZE= IN.STRINGLENGTH/2+50;}
	if (BP.PLUSSTRAND) BP.PLUSSTRAND=IN.STRINGLENGTH;
	IN.Positive_Head=IN.Tag_Copy;


/*	FILE* Original_File=File_Open(F.BINFILE,"rb");
	Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	if(!fread(Original_Text,Get_File_Size(Original_File),1,Original_File))fprintf (stderr,"Init(): Error loading genome...\n");
*/
	if (IN.STRINGLENGTH < SW_MIN_MATCH_LEN) SW_MIN_MATCH_LEN=IN.STRINGLENGTH-1;
	if (BP.MISHITS) 
	{
		Mishit_File1=File_Open(BP.MISFILE1,"w"); 
		Mishit_File2=File_Open(BP.MISFILE2,"w"); 
	}

	if(MISC_VERB)
	{
		if (R.COMPRESS) fprintf (stderr,"Using compressed index...\n");
		fprintf(stderr,"Single-end Alignment\n");
		if (BP.SMITH_WATERMAN)
		{
			fprintf(stderr,"Mate rescue mode : Flank size - %d\n",BP.SW_FLANKSIZE);
		}
		fprintf(stderr,"Max Gap : %d\n",BP.INSERTSIZE+BP.FLANKSIZE);
		//fprintf(stderr,"Using Pairing index %s\t %s\n", F.INDFILE,F.BLKFILE); 
		//fprintf (stderr,"------------------------------------------------------------------------------------------------------------------\n");
	}

}

void FreeQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t ) 
{
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> tmp; 
	while(!t.empty())
	{
		t.pop();
	}
	std::swap(t,tmp );
}

void Get_top1M_Alignments(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t_n,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & t, int & cutoff ) 
{
	if(t.size() < 1) 
	{
		return;
	}
	if(t.size()>1000)
		cutoff=90;
	else if(t.size()>500)
		cutoff=200;
	else if(t.size()>100)
		cutoff=300;
	else cutoff=500;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> tmp; 
	Alignment A=t.top();
	while(!t.empty() )
	{
	//if(labs(A.Loc-2853918318)<500) printf("\nllll %ld %d\n",A.Loc,A.Mismatch);
		if(labs(A.Score)>cutoff) break;
		if(A.Mismatch>MAX_MISMATCHES+1) 
		{
			t.pop();A=t.top();
			continue;
		}
		tmp.push(A);
		t.pop();
		A=t.top();
	}
	std::swap(t,tmp );
	t_n=t;
	FreeQ(tmp);
}

void FreeQ_Pair(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair>  & t ) 
{
	std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> tmp; 
	while(!t.empty())
	{
		t.pop();
	}
	std::swap(t,tmp );
}
void Get_Be_Alignments(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair>  & t ) 
{
	std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> tmp; 
	while(!t.empty())
	{
		t.pop();
	}
	std::swap(t,tmp );
}

const int FHSIZE=10;//00;
const int FHSKIP=500;

bool Unique_Hits(int Plus_Hits,int Minus_Hits,SARange & P,SARange & M)
{
	if(Plus_Hits)
	{
		if (P.Start!=P.End) 
		{
			return false; 
		}
		else 
		{
			return true;
		}
	}
	else
	{
		if (M.Start!=M.End) 
		{
			return false;
		}	
		else 
		{
			return true;
		}
	}
}

bool Check_Subopt(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,float & Top_BQScore,float & Sub_Score,float & Sub_BQScore, int & Quality_Score)
{
	if (SUPOPT_EXIST_FILTER) return false;
	else if (PHRED_FILTER) return Phred_Check(Plus_Hits,Minus_Hits,Top_Mis,Subopt_Mis,R,StringLength,MF,MC,Top_Score,Top_BQScore,Sub_Score,Sub_BQScore,Quality_Score);
	else SNP_Check(Plus_Hits,Minus_Hits,Top_Mis,Subopt_Mis,R,StringLength,MF,MC);
}

Alignment Realign(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength,READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments)
{
	s_align* Aln;
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	char Org_String[ORGSTRINGLENGTH],Cigar[MAX_SIGLEN];
	char Org_String_Ori[ORGSTRINGLENGTH];
	Cigar_Info Cig_Info;

	char* Current_Tag=(H.Sign=='+')? Read.Forward : Read.Complement;
	char* Current_Tag_raw=(H.Sign=='+')? Read.Forward_raw : Read.Complement_raw;
	A.Realigned=NO;
	A.Extend=0;
	int Jump=0;if(H.Sign=='-') Jump= 0+JUMP;
	s_profile* p = ssw_init((int8_t*)Current_Tag, StringLength, mata, n, 1);
	if((long)(H.Org_Loc+Jump)<0) H.Org_Loc = -Jump+1;
	if(H.hitType == '1')
		Get_Bases(Original_Text_CT,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String);
	else if(H.hitType == '4')
		Get_Bases(Original_Text_GA,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String);
	else {
		fprintf(stderr, "\nRealign Error!\n");
	}
	Get_Bases(Original_Text_Ori,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String_Ori);//momo

	Aln=mengyao_ssw_core(Org_String/*+Jump*/,StringLength, Current_Tag,StringLength+INDELGAP,0,0/*DP*/, p);
	if(Aln->score1 >= ACC_SCORE)
	{
		A.SW_Score=Aln->score1;
		ssw_cigar_processQ(H.Sign,Current_Tag_raw,Org_String_Ori,Aln,Cig_Info,Org_String,Aln->ref_begin1,Current_Tag,Aln->read_begin1,StringLength,R.Quality,NULL,0,0);
		if(Cig_Info.Mis > MAX_MISMATCHES  && false) 
		{
			//align_destroy(Aln);
			init_destroy(p); 
			A.Score==INT_MAX;
			A.Loc=INT_MAX;
			return A;
		}
		A.Loc=H.Org_Loc+Aln->ref_begin1;
		A.Score= -Cig_Info.Score;
		A.QScore=Cig_Info.QScore;
		A.BQScore=Cig_Info.BQScore;
		A.Mismatch=Cig_Info.Mis;
		A.Indel=Cig_Info.Indel_Count;
		A.swscore = Cig_Info.swscore;
		A.hitType = H.hitType;
		A.Clip=Aln->read_begin1;
		A.Sign=H.Sign;
		if(A.Indel)
		{
			if(!Dont_Push_To_Q)
			{
				Good_Alignments.push(A);
			}
		}
		else
		{
			//if(A.Mismatch> Inter_MM && !Dont_Push_To_Q)
			{
				Good_Alignments.push(A);
			}
		}
	}
	align_destroy(Aln);
	init_destroy(p); 
	return A;
	
}


//TODO
float Calc_Top_Score(MEMX & MF,MEMX & MC,float & Top_BQ,int Top_Mis,int StringLength,int Plus_Hits,int Minus_Hits,READ & R)
{
	if(Plus_Hits)
	{
		Top_BQ=BQSumX(MF.Hit_Array[0],R.Quality,Top_Mis,StringLength);
		return QSumX(MF.Hit_Array[0],R.Quality,Top_Mis,StringLength);
	}
	else
	{
		char Rev_Qual[100];
		Reverse_Quality(Rev_Qual,R,StringLength);
		Top_BQ=BQSumX(MC.Hit_Array[0],Rev_Qual,Top_Mis,StringLength);
		return QSumX(MC.Hit_Array[0],Rev_Qual,Top_Mis,StringLength);
	}
}

//------------------------
Alignment RealignFast(READ & RawR,char source,unsigned char* Original_Text,Hit_Info &  H,int StringLength, READ & R,int OFF,int Filter,bool Do_Filter)
{
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	if(!FASTSW)
	{
		A.Score=INT_MAX;
		return A; 
	}
	char Org_String[ORGSTRINGLENGTH];
	char Org_String_Ori[ORGSTRINGLENGTH];
	char *Quality,Rev_Qual[MAXTAG];
	Cigar_Info Cig_Info;


	char Real_String[R.Real_Len],*Current_Tag;
	char *Current_Tag_raw;
	//if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;RawR.Tag_Copy[RawR.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		if(H.Sign=='+')
		{
			Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
			Quality=R.Quality;
		}
		else
		{
			assert(false);
		}
		StringLength=R.Real_Len;
	}
	Current_Tag=R.Tag_Copy;
	Current_Tag_raw=RawR.Tag_Copy;


	Get_Bases(Original_Text_Ori,H.Org_Loc,StringLength+INDELGAP,Org_String_Ori);
	if(H.hitType == '1')
		Get_Bases(Original_Text_CT,H.Org_Loc,StringLength+INDELGAP,Org_String);
	else if(H.hitType == '4')
		Get_Bases(Original_Text_GA,H.Org_Loc,StringLength+INDELGAP,Org_String);
	else {
		fprintf(stderr, "\nRealign Fast Error!\n");
	}

	for(int i=0;i<StringLength+INDELGAP;i++)
	{
		*(Org_String+i)="ACGT"[*(Org_String+i)];
	}
	*(Org_String+StringLength+INDELGAP)=0;
	
	for(int i=0;i<StringLength+INDELGAP;i++)
	{
		*(Org_String_Ori+i)="ACGT"[*(Org_String_Ori+i)];
	}
	*(Org_String_Ori+StringLength+INDELGAP)=0;
	
	A=Fast_SW(Current_Tag_raw,source,Org_String,Current_Tag,Quality,OFF,'+',Org_String_Ori,R);
	if(A.SW_Score<match*9*StringLength/10) 
		A.Score=INT_MAX;
	A.Loc=H.Org_Loc+1;
	A.Sign=H.Sign; A.hitType = H.hitType;
	return A;
	
}

Alignment RealignX(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength, READ & R,bool Dont_Push_To_Q,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,char* Cigar,int & Clip_H,int & Clip_T,int & Filter,bool Do_Filter)
{
	s_align* Aln;
	Alignment A;
	A.Clip_H=A.Clip_T=INT_MAX;
	char Org_String[ORGSTRINGLENGTH];
	char Org_String_Ori[ORGSTRINGLENGTH];
	char *Quality,Rev_Qual[MAXTAG];
	Cigar_Info Cig_Info;

	char Real_String[R.Real_Len],*Current_Tag;
	char Real_String_raw[RawR.Real_Len],*Current_Tag_raw;
//printf("\n%s\n%s\n", R.Tag_Copy, RawR.Tag_Copy);
//	if(R.batmeth1) StringLength=R.Real_Len;
	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		RawR.Tag_Copy[RawR.Real_Len]=0;
		RawR.Quality[RawR.Real_Len]=0;
		if(H.Sign=='+')
		{
			Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
			Read2Bin(Real_String_raw,RawR.Tag_Copy,R.Real_Len);
			Quality=R.Quality;
			//H.Org_Loc-=JUMP-1; //Location changed by qw
		}
		else
		{
			Read2RevCBin(Real_String,R.Tag_Copy,R.Real_Len);
			Read2RevCBin(Real_String_raw,RawR.Tag_Copy,R.Real_Len);
			Reverse_Quality(Rev_Qual,R,R.Real_Len);
			//H.Org_Loc-=(R.Real_Len-StringLength)+INDELGAP-1;
			H.Org_Loc-=INDELGAP-1; //Location changed by qw
			//Offset=R.Real_Len-StringLength;
			Quality=Rev_Qual;
		}
		StringLength=R.Real_Len;
	}
	Current_Tag=Real_String;Current_Tag_raw=Real_String_raw;

//for(int i=0;i<StringLength;i++)
//	printf("%d",Current_Tag[i]);
	if(R.Real_Len < StringLength) StringLength = R.Real_Len;

	int Jump=0;if(H.Sign=='-') Jump=0 +JUMP;
	s_profile* p = ssw_init((int8_t*)Current_Tag, StringLength, mata, n, 1);//StringLength
	if( (long)(H.Org_Loc+Jump) <= 0) H.Org_Loc = (-Jump+1);
	Get_Bases(Original_Text_Ori,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String_Ori);//momo
	if(H.hitType == '1')
		Get_Bases(Original_Text_CT,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String);
	else if(H.hitType == '4')
		Get_Bases(Original_Text_GA,H.Org_Loc+Jump,StringLength+INDELGAP,Org_String);
	else {
		fprintf(stderr, "\nNo possible %c\n", H.hitType);
		exit(0);
	}

//
/*
printf("\ngenome %ld %d %d %d %c\n",H.Org_Loc,Jump,StringLength,INDELGAP,H.Sign);
for(int i=0;i<StringLength+INDELGAP;i++)
        printf("%d",Org_String_Ori[i]);
printf("\n");
for(int i=0;i<StringLength+INDELGAP;i++)
	printf("%d",Org_String[i]);
printf("\n====\n");
for(int i=0;i<StringLength;i++)
        printf("%d",Current_Tag[i]);
printf("\n");
for(int i=0;i<StringLength;i++)
        printf("%d",Current_Tag_raw[i]);
*/
//
	Aln=mengyao_ssw_core(Org_String,StringLength,Current_Tag,StringLength+INDELGAP,Filter,0/*DP*/, p);
//Ann_Info Ann; 
//Location_To_Genome(H.Org_Loc,Ann);
//printf("==== %s\n",Ann.Name);

//if( (labs(A.Loc-991289706)<=10 || labs(A.Loc-12218526)<=10) && !strcmp(Ann.Name,"chr19")) printf("\nAAAAAA %d\n",A.Loc);
	//if(Aln->score1 >= ACC_SCORE)
//printf("\nAAAAAA %d %d %s %d\n",Aln->mismatch_count, Aln->gap_count,A.Cigar, H.Org_Loc);
	if(Aln->score1 >= Filter)
	{
		A.Clip_H=Aln->read_begin1;A.Clip_T=0;
		if(Aln->read_end1!=StringLength-1) A.Clip_T=StringLength-1-Aln->read_end1;

		A.SW_Score=Aln->score1;
		ssw_cigar_processQ(H.Sign,Current_Tag_raw,Org_String_Ori,Aln,Cig_Info,Org_String,Aln->ref_begin1,Current_Tag,Aln->read_begin1,StringLength,Quality,A.Cigar,A.Clip_H,A.Clip_T);//,R
//printf("\nNNNNCigar %s %d %d %d %d %d %d %d %d %d",A.Cigar,Cig_Info.Score,Cig_Info.Mis,StringLength, Aln->read_begin1, Aln->ref_begin1, Aln->score1, Cig_Info.swscore, A.Clip_H, A.Clip_T);
		A.Loc=H.Org_Loc+Jump+Aln->ref_begin1;//+Offset;
		A.Score= -Cig_Info.Score;
		A.QScore=Cig_Info.QScore;
		A.BQScore=Cig_Info.BQScore;
		A.Mismatch=Cig_Info.Mis;
		A.Indel=Cig_Info.Indel_Count;
		A.Clip=A.Clip_H+A.Clip_T;
		A.swscore = Cig_Info.swscore;
		A.hitType = H.hitType;
		A.Sign=H.Sign;
		A.Realigned=1;
		A.Extend=1;
		if(Do_Filter && BOOST && Aln->score1 >Filter) 
		{
			Filter=Aln->score1;
		}
		if(!Dont_Push_To_Q)
		{
			Good_Alignments.push(A);
		}
		
	}
	else
	{
		A.Loc=INT_MAX;
		if(BOOST && Aln->score1>=0)
		{
			A.Score= -INT_MIN;
			A.Realigned=1;
			A.Extend=1;
			if(!Dont_Push_To_Q)
			{
				Good_Alignments.push(A);
			}
		}
	}
	align_destroy(Aln);
	init_destroy(p); 
	return A;
	
}

inline void Copy_A_to_H(Alignment & A,Hit_Info & H)
{
	H.Mismatch=A.Mismatch;
	H.Indel=A.Indel;
	H.Clip=A.Clip;
	H.Score= -A.Score;
	H.QScore= A.QScore;
	H.BQScore= A.BQScore;
	H.SW_Score= A.SW_Score;
	H.swscore = A.swscore;
	H.hitType = A.hitType;
	H.Loc=H.Org_Loc= A.Loc-1;
	H.Clip_H=A.Clip_H;H.Clip_T=A.Clip_T;
	if(A.Cigar[0])
		strcpy(H.Cigar,A.Cigar);
	else
	{
		H.Cigar[0]=0;		
		H.Cigar[1]='Z';
		A.Cigar[1]='Z';
	}

}

unsigned SA2Loc(BWT *revfmi,SARange S,int Pos,unsigned Conversion_Factor)
{
	if (S.Start==S.End) 
	{
		return S.Start;
	}
	else
	{
		assert (S.Start);
		return Conversion_Factor-BWTSaValue(revfmi,S.Start+Pos);
	}
	
}

void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T)
{
	Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
	int Thread_Num=0;
	pthread_attr_t Attrib;
	pthread_attr_init(&Attrib);
	pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);

	for (int i=0;i<NTHREAD;i++)
	{
		T.ThreadID=i;
		Thread_Info[i].Arg=T;
		//if(!(Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,NULL,Map_t,(void*) &Thread_Info[i].Arg))) Thread_Num++;
		Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib,Map_t,(void*) &Thread_Info[i].Arg);
		if(Thread_Info[i].r) {fprintf(stderr,"Launch_Threads():Cannot create thread..\n");exit(-1);} else Thread_Num++;
	}
	fprintf(stderr,"%d Threads runnning ...\n",Thread_Num);
	pthread_attr_destroy(&Attrib);

	for (int i=0;i<NTHREAD;i++)
	{
		pthread_join(Thread_Info[i].Thread,NULL);
	}
}

void Mode_Parameters(BATPARAMETERS BP)
{
	if (BP.SCANMODE==VERYSENSITIVE)
	{
		MODE=VERYSENSITIVE;
		CUT_MAX_SWALIGN=INT_MAX;
	}
	else if(BP.SCANMODE==SENSITIVE)
	{
		MODE=SENSITIVE;
		MAX_SW=20000;
		CUT_MAX_SWALIGN=20000;
	}
	else if(BP.SCANMODE==FAST)
	{
		MODE=FAST;
		CUT_MAX_SWALIGN=200;
	}
	else if(BP.SCANMODE==VERYFAST)
	{
		MODE=VERYFAST;
		//MAX_SW=200;
		CUT_MAX_SWALIGN=200;
	}

}

void Set_Force_Indel(bool & Force_Indel,int Last_Mis,Hit_Info & H,int Avg_Q)
{
	if(H.Indel || Last_Mis>1)
	{
		if(MODE<=FAST)
		{
			if(H.Indel+Last_Mis>1) 
			{
				if(Last_Mis)
				{
					int Avg_Score=(H.QScore)/(Last_Mis);
					if(Avg_Score>=Avg_Q)
						Force_Indel=true;
				}
			}
		}
		else
		{
			int Mis_Limit;
			if(MODE==VERYSENSITIVE) Mis_Limit=0; else Mis_Limit=1;
			if(Last_Mis>Mis_Limit) 
			{
				if(Last_Mis)
				{
					int Avg_Score=(H.QScore)/(Last_Mis);
					if(Avg_Score>=Avg_Q)
						Force_Indel=true;
				}
			}
			if(H.Indel)
				Force_Indel=true;
		}
	}
}

int Calculate_Average_Quality(READ & R)
{
	int Avg=0;
	for(int i=0;i<R.Real_Len;i++)
	{
		Avg+=R.Quality[i]-QUALITYCONVERSIONFACTOR;
	}
	Avg=Avg/R.Real_Len;
	return Avg/3;
}


void Set_Affinity()
{
	cpu_set_t Set;
	CPU_ZERO(&Set);

	if(sched_getaffinity(0,sizeof(cpu_set_t),&Set)<0)
	{
		fprintf(stderr,"Affinity could not be get..\n");
	}
	else
	{
		for (int i=0;i<CPU_SETSIZE;i++)
		{
			if(CPU_ISSET(i,&Set))
			{
				fprintf(stderr,"Bound to %d\n",i);
				CPU_ZERO(&Set);
				CPU_SET(i, &Set);
				if(sched_setaffinity(0, sizeof(Set), &Set)<0)
				{
					fprintf(stderr,"Affinity could not be set..\n");
				}
				return;
			}
		}
	}

}

bool Report_Single(READ & RawR,unsigned char* Original_Text,READ & R,Final_Hit & Single_File,const int StringLength,BATREAD & Read,bool & Print_Status,int Clip_H,int Clip_T,Alignment & A)
{
	Ann_Info Ann;
	Hit_Info H;
	assert(A.QScore!=INT_MAX);
	H.Org_Loc=A.Loc;H.Loc = A.Loc;H.Sign=A.Sign;H.QScore=A.QScore;H.Status=UNIQUEHIT;H.Indel=A.Indel;H.Clip=A.Clip;H.BQScore=A.BQScore;
	H.hitType = A.hitType;
	//A.Score= -A.Score;
	Location_To_Genome(H.Loc,Ann);H.Chr=Ann.Name;
	assert(A.Realigned);
	if (H.Loc+StringLength <= Ann.Size)
	{
		if(H.hitType == '1')
			Cigar_Check_And_Print(RawR,Original_Text_CT,H,Read,StringLength,Single_File,R,true,30,A,A.Clip_H,A.Clip_T,A.Cigar);
		else if(H.hitType == '4')
			Cigar_Check_And_Print(RawR,Original_Text_GA,H,Read,StringLength,Single_File,R,true,30,A,A.Clip_H,A.Clip_T,A.Cigar);
		else {
			fprintf(stderr, "\nUne error 5!\n"); exit(0);
		}
		Print_Status=true;
	}
}

void Alignments_RD(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A)
{
	std::map<unsigned,Alignment> D;
        while(!A.empty())
        {
                Alignment Aln=A.top();A.pop();
                std::map<unsigned,Alignment>::iterator I;
                if((I=D.find(Aln.Loc))==D.end())
                {
                        D[Aln.Loc]=Aln;
                }
                else if((I->second).Score < Aln.Score)
                {
                        D[Aln.Loc]=Aln;
                }
        }

        for(std::map<unsigned,Alignment>::iterator I=D.begin();I!=D.end();I++)
        {
		A.push(I->second);
	}
}
//Head 75bp
void Map_One_SEG_Head(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE)
{
	if(R.Real_Len<=75) return;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[SEG_SIZE],T_EOQ=R.Quality[SEG_SIZE];
	char T_EOS_raw=RawR.Tag_Copy[SEG_SIZE],T_EOQ_raw=RawR.Quality[SEG_SIZE];
	
	SEG_SIZE=75;SHIFT_SEG_SIZE=0;
	if(SEG_SIZE)
	{
		//for(int i=0;i<=SHIFT_SEG_SIZE;i++)
		for(int i=0;i<SEG_SIZE && (i+SHIFT_SEG_SIZE) < RawR.Real_Len;i++)
		{
			assert((R.Tag_Copy[i+SHIFT_SEG_SIZE]>='A' && R.Tag_Copy[i+SHIFT_SEG_SIZE]<='t'));
			R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_SEG_SIZE];R.Quality[i]=R.Quality[i+SHIFT_SEG_SIZE];
			RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_SEG_SIZE];RawR.Quality[i]=RawR.Quality[i+SHIFT_SEG_SIZE];
		}
		R.Tag_Copy[SEG_SIZE]='\0';
		RawR.Tag_Copy[SEG_SIZE]='\0';
	}
	
	R.Real_Len=0;
	for(;R.Tag_Copy[R.Real_Len]!='\0' && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
	RawR.Real_Len=R.Real_Len;
	
	//if(Segment_Length)
	//	R.Real_Len=Segment_Length;

	Conversion_Factor=revfmi->textLength-L.STRINGLENGTH;
	IN.Positive_Head=R.Tag_Copy;
	R.Read_Number=Actual_Tag;
	/*
	Process_Read(RawR,R,B,MF,MC);
	MF.Strand='+';MF.Larger_Than_Ten=0;MC.Strand='-';MC.Larger_Than_Ten=0;
	MF.Extend=MC.Extend=false;
	*/

	//
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
	
	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
	int Last_Mis;
	int Head_Top_Count;//# of top hits in H/T
	Quality_Score=QUAL_UNMAPPED;H.Status=UNMAPPED;
	
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
	int Paired_Score=INT_MAX;
	Align_Map_One(R,B,source,mismatch,L,fwfmi,revfmi,MF,MC,Paired_Score,RawR,Align_Hits,Alignments_tmp, 0, RTemp.Real_Len - R.Real_Len);
	Good_Alignments_tmp=Alignments_tmp;
	//Do_Mismatch_Scan(AM,'1',RawR,source,Original_Text,MF,MC,L,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Head_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
	int Loc_Mod = RTemp.Real_Len-SEG_SIZE; //-1;
	Fix_Offset_Head(Alignments_tmp,Alignments,Loc_Mod+JUMP,Loc_Mod+JUMP);
        Fix_Offset_Head(Good_Alignments_tmp,Good_Alignments,Loc_Mod+JUMP,Loc_Mod+JUMP);
	R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
    RawR=RawRTemp;
                
}
//Tail 75bp
void Map_One_SEG_Tail(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF2,MEMX & MC2,LEN & L,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int Segment_Length,int SEG_SIZE,int SHIFT_SEG_SIZE)
{
	if(R.Real_Len<=75) return;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[SEG_SIZE],T_EOQ=R.Quality[SEG_SIZE];
	char T_EOS_raw=RawR.Tag_Copy[SEG_SIZE],T_EOQ_raw=RawR.Quality[SEG_SIZE];
	
	//Tail75
	SEG_SIZE=75;SHIFT_SEG_SIZE=R.Real_Len-SEG_SIZE; //-1
	if(SEG_SIZE)
	{
    		//for(int i=0;i<=SHIFT_SEG_SIZE;i++)
		for(int i=0;i<SEG_SIZE && (i+SHIFT_SEG_SIZE) < R.Real_Len;i++)
		{
        		assert((R.Tag_Copy[i+SHIFT_SEG_SIZE]>='A' && R.Tag_Copy[i+SHIFT_SEG_SIZE]<='t'));
                	R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_SEG_SIZE];
			R.Quality[i]=R.Quality[i+SHIFT_SEG_SIZE];
             		RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_SEG_SIZE];
			RawR.Quality[i]=RawR.Quality[i+SHIFT_SEG_SIZE];
		}
	    	R.Tag_Copy[SEG_SIZE]='\0';
    		RawR.Tag_Copy[SEG_SIZE]='\0';
    	}

	R.Real_Len=0;
	for(;R.Tag_Copy[R.Real_Len]!='\0' && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
	RawR.Real_Len=R.Real_Len;
	Conversion_Factor=revfmi->textLength-L.STRINGLENGTH;

    if(Segment_Length)
	    R.Real_Len=Segment_Length;

/*		
    Process_Read(RawR,R,B,MF2,MC2);
    MF2.Strand='+';MF2.Larger_Than_Ten=0;MC2.Strand='-';MC2.Larger_Than_Ten=0;
    MF2.Extend=MC2.Extend=false;
*/    
    	//
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
	
	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
    int Tail_Top_Count;int Last_Mis;

	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
    	int Paired_Score=INT_MAX;
	Align_Map_One(R,B,source,mismatch,L,fwfmi,revfmi,MF2,MC2,Paired_Score,RawR,Align_Hits,Alignments_tmp, RTemp.Real_Len - R.Real_Len, 0);
	Good_Alignments_tmp=Alignments_tmp;
	//Do_Mismatch_Scan(AM,'2',RawR,source,Original_Text,MF2,MC2,L,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Tail_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
	Fix_Offset_Tail(Alignments_tmp,Alignments,SHIFT_SEG_SIZE+JUMP,SHIFT_SEG_SIZE+JUMP);
	Fix_Offset_Tail(Good_Alignments_tmp,Good_Alignments,SHIFT_SEG_SIZE+JUMP,SHIFT_SEG_SIZE+JUMP);

	R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
    RawR=RawRTemp;

}
//Head Half
void Map_One_Half_Head(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLH,MEMX & MCLH,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score)
{
	if(R.Real_Len<=SMALLSEED) return;
	int Half=SMALLSEED,SHIFT_Half=0;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[Half],T_EOQ=R.Quality[Half];
	char T_EOS_raw=RawR.Tag_Copy[Half],T_EOQ_raw=RawR.Quality[Half];
	
	//50 head
	for(int i=0;i<Half && (i+SHIFT_Half) < R.Real_Len;i++)
	{    
		assert((R.Tag_Copy[i+SHIFT_Half]>='A' && R.Tag_Copy[i+SHIFT_Half]<='t'));
        R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_Half];
		R.Quality[i]=R.Quality[i+SHIFT_Half];
        RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_Half];
		RawR.Quality[i]=RawR.Quality[i+SHIFT_Half];
	}
	R.Tag_Copy[Half]='\0';
	RawR.Tag_Copy[Half]='\0';
	R.Real_Len=RawR.Real_Len=Half;
	B.StringLength=Half;B.NCount=0;B.IGNOREHEAD=0;
/*
	Process_Read(RawR,R,B,MFLH,MCLH);
	MFLH.Strand='+';MFLH.Larger_Than_Ten=0;MCLH.Strand='-';MCLH.Larger_Than_Ten=0;
	MFLH.Extend=MCLH.Extend=false;
	MFLH.L=MCLH.L=L_Half;
*/	
		//
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
	Conversion_Factor=revfmi->textLength-L_Half.STRINGLENGTH;

	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
	
	int Head_Top_Count;int Last_Mis;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
	int Paired_Score=INT_MAX;
	Align_Map_One(R,B,source,mismatch,L_Half,fwfmi,revfmi,MFLH,MCLH,Paired_Score,RawR,Align_Hits,Alignments_tmp, 0, RTemp.Real_Len - R.Real_Len);
	Good_Alignments_tmp=Alignments_tmp;
	
	//Do_Mismatch_Scan(AM,'3',RawR,source,Original_Text,MFLH,MCLH,L_Half,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Head_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
	SHIFT_Half=RTemp.Real_Len-Half;
	//if(Alignments_tmp.size()>0) printf("\n%c %c %s\n",source,Alignments_tmp.top().Sign,R.Tag_Copy);
	Fix_Offset_Head(Alignments_tmp,Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
	Fix_Offset_Head(Good_Alignments_tmp,Good_Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);

	R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
    RawR=RawRTemp;
    
}
//Tail Half
void Map_One_Half_Tail(Align_Hit & Align_Hits,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLT,MEMX & MCLT,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score)
{
	if(R.Real_Len<=SMALLSEED) return;
	int Half=SMALLSEED,SHIFT_Half=R.Real_Len-Half;//-1;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[Half],T_EOQ=R.Quality[Half];
	char T_EOS_raw=RawR.Tag_Copy[Half],T_EOQ_raw=RawR.Quality[Half];
	
	//Tail 50
	for(int i=0;i<Half && (i+SHIFT_Half) < R.Real_Len;i++)
	{
		assert((R.Tag_Copy[i+SHIFT_Half]>='A' && R.Tag_Copy[i+SHIFT_Half]<='t'));
		R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_Half];
		R.Quality[i]=R.Quality[i+SHIFT_Half];
		RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_Half];
		RawR.Quality[i]=RawR.Quality[i+SHIFT_Half];
	}
	R.Tag_Copy[Half]='\0';
	RawR.Tag_Copy[Half]='\0';
	R.Real_Len=RawR.Real_Len=Half;
	B.StringLength=Half;B.NCount=0;B.IGNOREHEAD=0;
	/*
	Process_Read(RawR,R,B,MFLT,MCLT);
	MFLT.Strand='+';MFLT.Larger_Than_Ten=0;MCLT.Strand='-';MCLT.Larger_Than_Ten=0;
	MFLT.Extend=MCLT.Extend=false;
	MFLT.L=MCLT.L=L_Half;
	*/
			//
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
	Conversion_Factor=revfmi->textLength-L_Half.STRINGLENGTH;

	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
	int Tail_Top_Count;int Last_Mis;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
	int Paired_Score=INT_MAX;
	Align_Map_One(R,B,source,mismatch,L_Half,fwfmi,revfmi,MFLT,MCLT,Paired_Score,RawR,Align_Hits,Alignments_tmp, RTemp.Real_Len - R.Real_Len, 0);
	Good_Alignments_tmp=Alignments_tmp;
	//Do_Mismatch_Scan(AM,'4',RawR,source,Original_Text,MFLT,MCLT,L_Half,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Tail_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
//	if(Alignments_tmp.size()>0) printf("\n%c %c %s %d\n",source,Alignments_tmp.top().Sign,R.Tag_Copy,Alignments_tmp.top().Loc);
	Fix_Offset_Tail(Alignments_tmp,Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
	Fix_Offset_Tail(Good_Alignments_tmp,Good_Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
//	if(Alignments.size()>0) printf("\n=2=\n%c %c %s %d\n",source,Alignments.top().Sign,R.Tag_Copy,Alignments.top().Loc);
//if(Alignments.size()>0) printf("\n%s\n%s\n%d %d %ld",R.Tag_Copy,RTemp.Tag_Copy,SHIFT_Half,Alignments.size(),Alignments.top().Loc);
	R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
    RawR=RawRTemp;
	
}
//Head Twenty-Five
void Map_One_Half_Head(std::set<unsigned> & AM,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLH,MEMX & MCLH,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score)
{
	int Half=25,SHIFT_Half=0;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[Half],T_EOQ=R.Quality[Half];
	char T_EOS_raw=RawR.Tag_Copy[Half],T_EOQ_raw=RawR.Quality[Half];
	
	//50 head
	for(int i=0;i<Half && (i+SHIFT_Half) < R.Real_Len;i++)
	{    
		assert((R.Tag_Copy[i+SHIFT_Half]>='A' && R.Tag_Copy[i+SHIFT_Half]<='t'));
        	R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_Half];
		R.Quality[i]=R.Quality[i+SHIFT_Half];
        	RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_Half];
		RawR.Quality[i]=RawR.Quality[i+SHIFT_Half];
	}
	R.Tag_Copy[Half]='\0';
	RawR.Tag_Copy[Half]='\0';
	R.Real_Len=RawR.Real_Len=Half;
     
	B.StringLength=Half;B.NCount=0;B.IGNOREHEAD=0;
	Process_Read(RawR,R,B,MFLH,MCLH);
	MFLH.Strand='+';MFLH.Larger_Than_Ten=0;MCLH.Strand='-';MCLH.Larger_Than_Ten=0;
	MFLH.Extend=MCLH.Extend=false;
	MFLH.L=MCLH.L=L_Half;
	
		//
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
//	Conversion_Factor=revfmi->textLength-L_Half.STRINGLENGTH;

	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
	
	int Head_Top_Count;int Last_Mis;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
//	Do_Mismatch_Scan(AM,'3',RawR,source,Original_Text,MFLH,MCLH,L_Half,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Head_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
	SHIFT_Half=RTemp.Real_Len-Half;
	//if(Alignments_tmp.size()>0) printf("\n%c %c %s\n",source,Alignments_tmp.top().Sign,R.Tag_Copy);
	Fix_Offset_Head(Alignments_tmp,Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
	Fix_Offset_Head(Good_Alignments_tmp,Good_Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);

	R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
    RawR=RawRTemp;
    
}
//Tail Thrity
void Map_One_Half_Tail(std::set<unsigned> & AM,int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MFLT,MEMX & MCLT,LEN & L_Half,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score)
{
	int Half=50,SHIFT_Half=R.Real_Len-Half;//-1;
	READ RTemp=R;BATREAD BTemp=B;
	READ RawRTemp=RawR;
	char T_EOS=R.Tag_Copy[Half],T_EOQ=R.Quality[Half];
	char T_EOS_raw=RawR.Tag_Copy[Half],T_EOQ_raw=RawR.Quality[Half];
	
	//Tail 50
	for(int i=0;i<Half && (i+SHIFT_Half) < R.Real_Len;i++)
	{
		assert((R.Tag_Copy[i+SHIFT_Half]>='A' && R.Tag_Copy[i+SHIFT_Half]<='t'));
		R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_Half];
		R.Quality[i]=R.Quality[i+SHIFT_Half];
		RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_Half];
		RawR.Quality[i]=RawR.Quality[i+SHIFT_Half];
	}
	R.Tag_Copy[Half]='\0';
	RawR.Tag_Copy[Half]='\0';
	R.Real_Len=RawR.Real_Len=Half;
	B.StringLength=Half;B.NCount=0;B.IGNOREHEAD=0;
	Process_Read(RawR,R,B,MFLT,MCLT);
	MFLT.Strand='+';MFLT.Larger_Than_Ten=0;MCLT.Strand='-';MCLT.Larger_Than_Ten=0;
	MFLT.Extend=MCLT.Extend=false;
	MFLT.L=MCLT.L=L_Half;
	
	FreeQ(Alignments);
	FreeQ(Good_Alignments);
	Conversion_Factor=revfmi->textLength-L_Half.STRINGLENGTH;

	if (R.NCount > NCOUNT) 
	{
		R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
		R=RTemp;B=BTemp;
		RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
		RawR=RawRTemp;
		return;
	}
	int Tail_Top_Count;int Last_Mis;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Good_Alignments_tmp;
//	Do_Mismatch_Scan(AM,'4',RawR,source,Original_Text,MFLT,MCLT,L_Half,fwfmi,revfmi,mismatch,mismatch,Last_Mis,Tail_Top_Count,H,Quality_Score,R,B,Mishit_File,Conversion_Factor,Alignments_tmp,Good_Alignments_tmp);
//	if(Alignments_tmp.size()>0) printf("\n%c %c %s %d\n",source,Alignments_tmp.top().Sign,R.Tag_Copy,Alignments_tmp.top().Loc);
	Fix_Offset_Tail(Alignments_tmp,Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
	Fix_Offset_Tail(Good_Alignments_tmp,Good_Alignments,SHIFT_Half+JUMP,SHIFT_Half+JUMP);
//	if(Alignments.size()>0) printf("\n=2=\n%c %c %s %d\n",source,Alignments.top().Sign,R.Tag_Copy,Alignments.top().Loc);
//if(Alignments.size()>0) printf("\n%s\n%s\n%d %d %ld",R.Tag_Copy,RTemp.Tag_Copy,SHIFT_Half,Alignments.size(),Alignments.top().Loc);
	R.Tag_Copy[Half]=T_EOS;R.Quality[Half]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[Half]=T_EOS_raw;RawR.Quality[Half]=T_EOQ_raw;
    RawR=RawRTemp;
	
}

///Quarter
void Map_One_Quart_Head(int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF,MEMX & MC,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Third,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int SEG_SIZE,int SHIFT_SEG_SIZE,Pair* & Pairs)
{
	if(R.Real_Len<=75) return;
        READ RTemp=R;BATREAD BTemp=B;
        READ RawRTemp=RawR;
        char T_EOS=R.Tag_Copy[SEG_SIZE],T_EOQ=R.Quality[SEG_SIZE];
        char T_EOS_raw=RawR.Tag_Copy[SEG_SIZE],T_EOQ_raw=RawR.Quality[SEG_SIZE];

        SEG_SIZE=75;SHIFT_SEG_SIZE=0;
        if(SEG_SIZE)
        {
                //for(int i=0;i<=SHIFT_SEG_SIZE;i++)
                for(int i=0;i<=SEG_SIZE && (i+SHIFT_SEG_SIZE) < RawR.Real_Len;i++)
                {
                        assert((R.Tag_Copy[i+SHIFT_SEG_SIZE]>='A' && R.Tag_Copy[i+SHIFT_SEG_SIZE]<='t'));
                        R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_SEG_SIZE];R.Quality[i]=R.Quality[i+SHIFT_SEG_SIZE];
                        RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_SEG_SIZE];RawR.Quality[i]=RawR.Quality[i+SHIFT_SEG_SIZE];
                }
                R.Tag_Copy[SEG_SIZE]='\0';
                RawR.Tag_Copy[SEG_SIZE]='\0';
        }

        R.Real_Len=0;
        for(;R.Tag_Copy[R.Real_Len]!='\0' && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
        RawR.Real_Len=R.Real_Len;

        //if(Segment_Length)
        //      R.Real_Len=Segment_Length;

        Conversion_Factor=revfmi->textLength-L.STRINGLENGTH;
        IN.Positive_Head=R.Tag_Copy;
        R.Read_Number=Actual_Tag;
        Process_Read(RawR,R,B,MF,MC);
        MF.Strand='+';MF.Larger_Than_Ten=0;MC.Strand='-';MC.Larger_Than_Ten=0;
        MF.Extend=MC.Extend=false;

        //

                READ R_Head;BATREAD B_Head;
                READ RawR_Head;
	int SEG=25;
	int Shift=75-SEG;
                for (int i=0;i<SEG;i++)
                {
                        R_Head.Tag_Copy[i]=R.Tag_Copy[i];
                        RawR_Head.Tag_Copy[i]=RawR.Tag_Copy[i];
                }
		RawR_Head.Tag_Copy[SEG]=0;R_Head.Tag_Copy[SEG]=0;
                B_Head.StringLength=SEG;B_Head.NCount=0;B_Head.IGNOREHEAD=0;
                Process_Read(RawR_Head,R_Head,B_Head,MFH,MCH);
                MFH.Strand='+';MFH.Larger_Than_Ten=0;MCH.Strand='-';MCH.Larger_Than_Ten=0;MFH.Extend=MCH.Extend=false;
                MFH.L=MCH.L=L_Third;

                READ R_Tail;BATREAD B_Tail;
                READ RawR_Tail;
                for (int i=0;i<SEG;i++)
                {
                        R_Tail.Tag_Copy[i]=R.Tag_Copy[i+Shift];
                        RawR_Tail.Tag_Copy[i]=RawR.Tag_Copy[i+Shift];
                }
		RawR_Tail.Tag_Copy[SEG]=R_Tail.Tag_Copy[SEG]=0;
                B_Tail.StringLength=SEG;B_Tail.NCount=0;B_Tail.IGNOREHEAD=0;
                Process_Read(RawR_Tail,R_Tail,B_Tail,MFT,MCT);
                MFT.Strand='+';MFT.Larger_Than_Ten=0;MCT.Strand='-';MCT.Larger_Than_Ten=0;MFT.Extend=MCT.Extend=false;
                MFT.L=MCT.L=L_Third;

	FreeQ(Alignments);
        FreeQ(Good_Alignments);

	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
	Do_Indel(source,RawR,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,MF,MC,MF,MC,MFH,MCH,MFT,MCT,L.STRINGLENGTH,Pairs,R,Alignments_tmp,H);

	int SHIFT_SIZE=R.Real_Len-75;
	Fix_Offset_Head(Alignments_tmp,Alignments,SHIFT_SIZE+JUMP,SHIFT_SIZE+JUMP);
	Good_Alignments=Alignments;

        R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
        R=RTemp;B=BTemp;
    RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
    RawR=RawRTemp;
}

void Map_One_Quart_Tail(int mismatch,READ & RawR,char source,RQINDEX & RQHALF,RQINDEX & RQ,unsigned char* Original_Text,unsigned Entries,BWT* fwfmi,BWT* revfmi,READ & R,BATREAD & B,unsigned & Conversion_Factor,MEMX & MF2,MEMX & MC2,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,LEN & L,LEN & L_Third,unsigned & Actual_Tag,Final_Hit &  Single_File,FILE* Mishit_File,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,bool PRINT,Hit_Info & H,int & Quality_Score,int SEG_SIZE,int SHIFT_SEG_SIZE,Pair* & Pairs)
{
	if(R.Real_Len<=75) return;
        READ RTemp=R;BATREAD BTemp=B;
        READ RawRTemp=RawR;
        char T_EOS=R.Tag_Copy[SEG_SIZE],T_EOQ=R.Quality[SEG_SIZE];
        char T_EOS_raw=RawR.Tag_Copy[SEG_SIZE],T_EOQ_raw=RawR.Quality[SEG_SIZE];

        //Tail75
        SEG_SIZE=75;SHIFT_SEG_SIZE=R.Real_Len-SEG_SIZE;//-1;
    if(SEG_SIZE)
    {
        //for(int i=0;i<=SHIFT_SEG_SIZE;i++)
        for(int i=0;i<SEG_SIZE && (i+SHIFT_SEG_SIZE) < R.Real_Len;i++)
        {
                assert((R.Tag_Copy[i+SHIFT_SEG_SIZE]>='A' && R.Tag_Copy[i+SHIFT_SEG_SIZE]<='t'));
            R.Tag_Copy[i]=R.Tag_Copy[i+SHIFT_SEG_SIZE];
                        R.Quality[i]=R.Quality[i+SHIFT_SEG_SIZE];
            RawR.Tag_Copy[i]=RawR.Tag_Copy[i+SHIFT_SEG_SIZE];
                        RawR.Quality[i]=RawR.Quality[i+SHIFT_SEG_SIZE];
                }
            R.Tag_Copy[SEG_SIZE]='\0';
        RawR.Tag_Copy[SEG_SIZE]='\0';
    }

        R.Real_Len=0;
        for(;R.Tag_Copy[R.Real_Len]!='\0' && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
        RawR.Real_Len=R.Real_Len;

//    if(Segment_Length)
  //          R.Real_Len=Segment_Length;

    Process_Read(RawR,R,B,MF2,MC2);
    MF2.Strand='+';MF2.Larger_Than_Ten=0;MC2.Strand='-';MC2.Larger_Than_Ten=0;
    MF2.Extend=MC2.Extend=false;

        //
	Conversion_Factor=revfmi->textLength-L.STRINGLENGTH;

                READ R_Head;BATREAD B_Head;
                READ RawR_Head;
        int SEG=25;
        int Shift=75-SEG;
                for (int i=0;i<SEG;i++)
                {
                        R_Head.Tag_Copy[i]=R.Tag_Copy[i];
                        RawR_Head.Tag_Copy[i]=RawR.Tag_Copy[i];
                }
		RawR_Head.Tag_Copy[SEG]=R_Head.Tag_Copy[SEG]='\0';
                B_Head.StringLength=SEG;B_Head.NCount=0;B_Head.IGNOREHEAD=0;
                Process_Read(RawR_Head,R_Head,B_Head,MFH,MCH);
                MFH.Strand='+';MFH.Larger_Than_Ten=0;MCH.Strand='-';MCH.Larger_Than_Ten=0;MFH.Extend=MCH.Extend=false;
                MFH.L=MCH.L=L_Third;

                READ R_Tail;BATREAD B_Tail;
                READ RawR_Tail;
                for (int i=0;i<SEG;i++)
                {
                        R_Tail.Tag_Copy[i]=R.Tag_Copy[i+Shift];
                        RawR_Tail.Tag_Copy[i]=RawR.Tag_Copy[i+Shift];
                }
		RawR_Tail.Tag_Copy[SEG]=R_Tail.Tag_Copy[SEG]='\0';
                B_Tail.StringLength=SEG;B_Tail.NCount=0;B_Tail.IGNOREHEAD=0;
                Process_Read(RawR_Tail,R_Tail,B_Tail,MFT,MCT);
                MFT.Strand='+';MFT.Larger_Than_Ten=0;MCT.Strand='-';MCT.Larger_Than_Ten=0;MFT.Extend=MCT.Extend=false;
                MFT.L=MCT.L=L_Third;

        FreeQ(Alignments);
        FreeQ(Good_Alignments);

        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> Alignments_tmp;
        Do_Indel(source,RawR,RQHALF,RQ,Original_Text,Entries,fwfmi,revfmi,MF2,MC2,MF2,MC2,MFH,MCH,MFT,MCT,L.STRINGLENGTH,Pairs,R,Alignments_tmp,H);

        int SHIFT_SIZE=R.Real_Len-75;
        Fix_Offset_Tail(Alignments_tmp,Alignments,SHIFT_SIZE+JUMP,SHIFT_SIZE+JUMP);
        Good_Alignments=Alignments;

        R.Tag_Copy[SEG_SIZE]=T_EOS;R.Quality[SEG_SIZE]=T_EOQ;
	R=RTemp;B=BTemp;
    RawR.Tag_Copy[SEG_SIZE]=T_EOS_raw;RawR.Quality[SEG_SIZE]=T_EOQ_raw;
    RawR=RawRTemp;
}


void Get_Basic_MapQ(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,Alignment & C1, int & MapQ)
{
	MapQ=1;
	if(!Good_Alignments.empty())
	{
		C1=Good_Alignments.top();C1.Score= -C1.Score;
		Good_Alignments.pop();
		if(!Good_Alignments.empty())
		{
			Alignment C2=Good_Alignments.top();C2.Score= -C2.Score;
			if(labs(C1.Score-C2.Score)<=10)
				MapQ=0;
		}
		//else
		//	MapQ=labs(C1.Score-C2.Score)+1;
		Good_Alignments.push(C1);
	}
	else
		MapQ= -1;
}

void Adjust_Alignments(READ & Raw,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp)
{
	if(Raw.Real_Len<75) return;
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> A;
	while(!Alignments.empty())
	{
		Pop_And_Realign(Raw,Original_Text,Alignments,A,Offset,RTemp,BTemp);
	}
	//int cutoff=(MAX_MISMATCHES+5)*30+10;
	//Get_top1M_Alignments(Alignments,A,cutoff);
	Alignments=A;

}

void Pop_And_Realign(READ & RawR,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,int Offset,READ & RTemp, BATREAD & BTemp)
{
	int Filter=0,ClipT,ClipH;
	char CIG2[MAX_SIGLEN];
	Hit_Info H;

	Alignment Aln=Alignments.top();
	Alignments.pop();
	if(Aln.Realigned==NO)
		Aln.Clip_H=Aln.Clip_T=Aln.Clip=0;
	assert(Aln.Realigned==NO ||Aln.Realigned==1);

//printf("\nExtend %d %s\n",Aln.Extend,Aln.Cigar);
	if(Aln.Extend==1)
	{//printf("\nKKKKKKK %s %d\n",Aln.Cigar,Aln.Score);
		A.push(Aln);
		return;
	}
	Aln.Extend=1;

	H.Org_Loc=Aln.Loc;H.Sign=Aln.Sign;
	if(H.Sign=='+')
	{
		H.Org_Loc-=Offset;
		if(Aln.Clip_H<=2)
			H.Org_Loc-=Aln.Clip_H;
	}
	else
	{
		H.Org_Loc+=Offset;
		if(Aln.Clip_T<=2)
			H.Org_Loc+=Aln.Clip_T;
	}
	H.hitType = Aln.hitType;
	RealignX(RawR,Original_Text,H,BTemp,75,RTemp,false,A,CIG2,ClipH,ClipT);
}

bool SW_List(READ & RawR,unsigned char* Original_Text,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Offset,READ & RTemp, BATREAD & BTemp,bool & List_Exceeded,int & LS)
{
	int Enum_Hits=0,Top_Hits=0;
	int List_Size=Good_Alignments.size();
	int Top_Score;
	LS+=List_Size;
	bool Top_Scanned=false;

	for(int i=0;Enum_Hits<MAX_PER_LIST && !Good_Alignments.empty();i++)
	{
		if(List_Size>=MAX_PER_LIST && !Top_Scanned)
		{
			Enum_Hits++;
			if(i==0)
				Top_Score=Good_Alignments.top().Score;
			else if(labs(Good_Alignments.top().Score-Top_Score)>10)
			{
				Top_Scanned=true;
				Enum_Hits=Top_Hits;
			}
				//break;
		}
		else
			Enum_Hits++;
		LS--;
		Pop_And_Realign(RawR,Original_Text,Good_Alignments,Alignments,Offset,RTemp,BTemp);
	}
	if(Enum_Hits==MAX_PER_LIST)
	{
		List_Exceeded=true;
		return true;
	}
	else
	{
		List_Exceeded=false;
		return false;
	}
}

bool Correct_Orientation(Alignment A,Alignment A_P,int Extra_Bit)
{
	if(A.Sign==A_P.Sign)
		return false;
	if(A.Sign=='+')
	{
		if(A_P.Sign!='-') return false;
		//assert(A_P.Sign=='-');
		//if(A.Loc<=A_P.Loc)
		{
			if(ESTIMATE)
				return true;
			if(labs(A_P.Loc-A.Loc)<=INSERTSIZE+3*STD+Extra_Bit || labs(A_P.Loc-A.Loc)<=550)
				return true;
		}
	}
	else
	{
		if( A_P.Sign!='+' ) return false;
		//assert(A_P.Sign=='+');
		//if(A.Loc>=A_P.Loc)
		{
			if(ESTIMATE)
				return true;
			if(labs(A.Loc-A_P.Loc)<=INSERTSIZE+3*STD+Extra_Bit || labs(A.Loc-A_P.Loc)<=550 ) //2 to 3 moxian
				return true;
		}
	}
	return false;
}

const int MAX_PROPER_PAIRS=2000;//moxian
bool Find_Paired(int & paired_score,bool & Unique,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & B,std::map<unsigned,Alignment> & D,std::map<unsigned,Alignment> & D_P,int Extra_Bit,int SW_Compare)
{
	Alignment Pairings[MAX_PROPER_PAIRS+1];int Pairings_Index=0;
	int Count1=0,Count2=0;
	while(!A.empty())
	{
		Count1++;
		Alignment Aln=A.top();A.pop();
		if(Aln.Loc==UINT_MAX) continue;
		std::map<unsigned,Alignment>::iterator I;
		if((I=D.find(Aln.Loc))==D.end())
		{
			Count2++;
			D[Aln.Loc]=Aln;
		}
		else if((I->second).Score < Aln.Score)
		{
			Count2++;
			D[Aln.Loc]=Aln;
		}
	}
	Count1=0,Count2=0;
	while(!B.empty())
	{
		Count1++;
		Alignment Aln=B.top();B.pop();
		std::map<unsigned,Alignment>::iterator I;
		if(Aln.Loc==UINT_MAX) continue;
		if((I=D_P.find(Aln.Loc))==D_P.end())
		{
			Count2++;
			D_P[Aln.Loc]=Aln;
		}
		else if((I->second).Score < Aln.Score)
		{
			Count2++;
			D_P[Aln.Loc]=Aln;
		}
	}

	Alignment Head,Tail;	
	Alignment Sub_Opt_Head,Sub_Opt_Tail;	
	int Sub_Opt_Score=INT_MAX;
	int Max_H_Score=INT_MIN,Max_T_Score=INT_MIN;
//	bool Unique=true;//unique best pair..
//printf("\n=-=-=-= %d %d\n",D.size(),D_P.size());
	for(std::map<unsigned,Alignment>::iterator I=D.begin();I!=D.end() && Pairings_Index<MAX_PROPER_PAIRS;I++)
	{
		std::map<unsigned,Alignment>::iterator Nearest_Pair=D_P.lower_bound(I->first-(INSERTSIZE+3*STD+Extra_Bit));

		while(Nearest_Pair!=D_P.end() && (labs(Nearest_Pair->first-I->first) < INSERTSIZE+3*STD+Extra_Bit))
		{//printf("\n%ld %ld %d %d %d\n",Nearest_Pair->first,I->first,Nearest_Pair->first-I->first,(Nearest_Pair->second).Score,(I->second).Score);
			int Paired_Score=(I->second).Score+(Nearest_Pair->second).Score;
			if(Max_H_Score<(I->second).Score)
			{
				Max_H_Score=(I->second).Score;
			}
			if(Max_T_Score<(Nearest_Pair->second).Score)
			{
				Max_T_Score=(Nearest_Pair->second).Score;
			}
			
			if (Correct_Orientation(I->second,Nearest_Pair->second,Extra_Bit))
			{
				if(!Pairings_Index)
				{
					Head=I->second;Tail=Nearest_Pair->second;
					Max_H_Score=Head.Score;Max_T_Score=Tail.Score;
//printf("\nCCCCCCC %ld %ld %d %d %s %s %d %d %d\n",Head.Loc,Tail.Loc,Head.Score,Tail.Score,Head.Cigar,Tail.Cigar,Head.Mismatch,Tail.Mismatch,Unique);
				}
				else if(Paired_Score >= (Head.Score+Tail.Score))
				{
					if(Paired_Score <= (Head.Score+Tail.Score)+10)
					{
						Unique=false;
					}
					else
					{
						Unique=true;
					}
					if(!Unique && (labs(Head.Loc-(I->second).Loc)<=20 || labs(Tail.Loc-(Nearest_Pair->second).Loc)<=20) ) Unique=true;
					Sub_Opt_Score=Head.Score+Tail.Score;
					Sub_Opt_Head=Head;Sub_Opt_Tail=Tail;
					Head=I->second;Tail=Nearest_Pair->second;
					paired_score=Head.Indel+Head.Mismatch+Tail.Indel+Tail.Mismatch;
//					printf("\nCCCCCCC %ld %ld %d %d %s %s %d %d %d\n",Head.Loc,Tail.Loc,Head.Score,Tail.Score,Head.Cigar,Tail.Cigar,Head.Mismatch,Tail.Mismatch,Unique);
				}
				Pairings_Index++;
			}
			else
			{
				if(Pairings_Index && Sub_Opt_Score==INT_MAX)//not the first hit..
				{
					Sub_Opt_Score=Head.Score+Tail.Score;
					Sub_Opt_Head=Head;Sub_Opt_Tail=Tail;
				}
				else if(Paired_Score>Sub_Opt_Score)
				{
					Sub_Opt_Score=Head.Score+Tail.Score;
					Sub_Opt_Head=Head;Sub_Opt_Tail=Tail;
				}
			}
			if(Head.Score+Tail.Score > (I->second).Score) break;
			if(Head.Score+Tail.Score > (Nearest_Pair->second).Score) D_P.erase(Nearest_Pair++);
			else Nearest_Pair++;
		}
		
//		if(Head.Score+Tail.Score > (I->second).Score) //because here is minus
//		{   printf("\nScore %d %d %d\n",Head.Score,Tail.Score,(I->second).Score);
//			D.erase(I++);
//		}else I++;
	}
//printf("\nHead %d %d %d %d\n",Head.Score,Head.Mismatch,Tail.Score,Tail.Mismatch);

	if(Sub_Opt_Score!=INT_MAX)
	{
		if(Head.Score<Max_H_Score)//Sub_Opt_Head.Score)
		{
			if(Unique)
				if(Head.Score < Sub_Opt_Head.Score) Sub_Opt_Head.Score=Head.Score-30;
			else
			{
				Sub_Opt_Tail.Score=Tail.Score;
				Sub_Opt_Head.Score=Head.Score;
			}
		}
		if(Tail.Score<Max_T_Score)//Sub_Opt_Tail.Score)
		{
			if(Unique)
				if(Tail.Score < Sub_Opt_Tail.Score) Sub_Opt_Tail.Score=Tail.Score-30;
			else
			{
				Sub_Opt_Head.Score=Head.Score;
				Sub_Opt_Tail.Score=Tail.Score;
			}
		}
	}
//printf("\nSub_Opt_Head %d %d %d %d\n",Sub_Opt_Head.Score,Sub_Opt_Head.Mismatch,Sub_Opt_Tail.Score,Sub_Opt_Tail.Mismatch);
	if(Pairings_Index)
	{
		if(false/*moxian*/ && Pairings_Index>=MAX_PROPER_PAIRS)
		{
			A.push(Head);B.push(Tail);
			A.push(Head);B.push(Tail);
			return true;
			//printf("MAX_PROPER_PAIRS limit exceeded..\n");
		}
		if(Unique)
		{
			A.push(Head);B.push(Tail);
			/*if(Sub_Opt_Score!=INT_MAX)
			{
				A.push(Sub_Opt_Head);B.push(Sub_Opt_Tail);
			}*///qw
		}
		else
		{
			assert(Sub_Opt_Score!=INT_MAX);
			A.push(Head);B.push(Tail);
			//A.push(Sub_Opt_Head);B.push(Sub_Opt_Tail); //qw
		}
	}
	else
	{
		return false;
	}
	return true;
}

bool Align_Difference(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,unsigned U)
{
	Alignment Aln=Alignments.top();//check if actual best alignment differes from the initial..
	if(ulabs(Aln.Loc,U)>75)
		return true;
	else
		return false;

}

unsigned ulabs(unsigned A,unsigned B)
{
	if(A>B)
	{
		return (A-B);
	}
	else
		return (B-A);
}

void Get_Head_Map(std::map<unsigned,Alignment> & D)
{
// get mismatch+Indel <= 1
	Alignment A;
	for(std::map<unsigned,Alignment>::iterator I=D.begin();I!=D.end();/*I++*/)
	{
		A=I->second;
		if(A.Mismatch+A.Indel > 1) //because here is minus
                {    
                        D.erase(I++);
                }else I++;
	}
}

void Mate_Rescue_New(int mismatch,int Indel,bool & Do_Rescue,bool & realign,std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,READ & RawR,READ & RawM,unsigned char* Original_Text,char source1,char source2,READ & RTemp,READ & RTemp_P,BATREAD & BTemp,BATREAD & BTemp_P,int Read_Length,int Read_Length2,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_tmp,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments_tmp_P,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Good_Alignments_P,Hit_Info & H1,Hit_Info & H1_P,FILE* Single_File,int Quality_Score1,int Quality_Score1_P,Alignment & A1,int MapQ2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi)
{
	assert(!ESTIMATE);
	if(Alignments_tmp.empty())
	{
		assert(false);
	}
	std::map<unsigned,Alignment> D,D_P;

	BTemp.StringLength=Read_Length;
	RTemp.Real_Len=Read_Length;
	if(source1!='0')
		Process_Read_Basic(RawR,RTemp,BTemp);

	Final_Hit Head_Hit,Mate_Hit;
	if(realign)
		Adjust_Alignments(RawR,Original_Text,Alignments_tmp,0,RTemp,BTemp);
	ALIGNMENT_Q Alignments=Alignments_tmp;
//printf("\n-=-=-= %s %ld %d %d || %d %d\n",Alignments_tmp.top().Cigar,Alignments_tmp.top().Loc,realign,RTemp.Real_Len,mismatch,Indel);

	{
		Remove_Dup_Top(Alignments,Read_Length);
		H1.Status=UNMAPPED;
		if(!Report_SW_Hits(RawR,source1,Original_Text,0,RTemp,Head_Hit,Read_Length,BTemp,H1,Quality_Score1,Alignments,Good_Alignments,0/*Force_Indel*/,true,false))
		{
			return;
		}
		Mate_Hit.Loc=INT_MAX;
		Print_Pair(Alignments_Reslut,H1,H1_P,H1.Sign,H1_P.Sign,RawR,RawM,Original_Text,source1,source2,Single_File,Head_Hit,Mate_Hit,RTemp,RTemp_P,MF,MC,fwfmi,revfmi);

		return;
	}

}

void Remove_Dup_Top(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  & Alignments,int Gap)
{
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment>  Alignnew;
	if(Alignments.size()<=1)
	{
		return;
	}
	assert(!Alignments.empty());

	Alignment Top_Aln=Alignments.top();
	Alignnew.push(Top_Aln);
	Alignments.pop();
	Alignment Sub_Aln=Alignments.top();

        while(1) {
                if(Alignments.empty()) break;
                if(!Alignments.empty() && ( dist(Sub_Aln.Loc,Top_Aln.Loc) < Gap)  )
                {
                        Alignments.pop();
                        Sub_Aln=Alignments.top();
                }else if(Alignments.size()>1){
                        Alignnew.push(Sub_Aln);
                        Top_Aln = Sub_Aln;
                        Alignments.pop();
                        Sub_Aln=Alignments.top();
                }else {
                        Alignnew.push(Alignments.top());
                        Alignments.pop();
                }
        }
        Alignments = Alignnew;
}

void get_len_match(Alignment & A, int & headM, int & TailM){
	const char* cig=A.Cigar;
	char temp[8];
	int numberS = 0, headS = 0;
	unsigned n=0;

	while(*cig!='\0')
	{
		if(*cig>='0' && *cig<='9')
		{
			temp[n]=*cig;
			cig++;n++;
		}else if(*cig=='S')
		{
			temp[n]='\0';int length=atoi(temp);
			if(numberS == 0) headS = length;
			numberS++;
			cig++;n=0;
		}else if(*cig=='M')
		{
			temp[n]='\0';int length=atoi(temp);
			if(numberS == 0 || (numberS ==1 && headS < 10) ) headM += length;
			else TailM += length;
			cig++;n=0;
		}else if(*cig=='I')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			cig++;n=0;
		}else if(*cig=='D')
		{
			int i;temp[n]='\0';int length=atoi(temp);
			cig++;n=0;
		}else
		{
			return;
		}
	}
}
bool Align_comp(Alignment & A, Alignment & B, int Aheadlen, int Ataillen, int Bheadlen, int Btaillen, int Real_len)
{
	if(A.Sign != B.Sign) {
		if((Ataillen > 20 && Aheadlen > 20 && Bheadlen > 20 && Btaillen>20 && dist(Ataillen, Btaillen)<=10) ||
		( Ataillen > 20 && Aheadlen > 20 && Bheadlen > 20 && Btaillen>20 && dist(Aheadlen, Bheadlen)<=10) ) {
			return true;
		}
	}else{
		if(Ataillen > 20 && Aheadlen > 20 && Bheadlen > 20 && Btaillen>20 && dist(Ataillen + Bheadlen, Real_len)<=10){
			return true;
		}
		if(Ataillen > 20 && Aheadlen > 20 && Bheadlen > 20 && Btaillen>20 && dist(Btaillen + Aheadlen, Real_len)<=10){
			return true;
		}
	}
	return false;
}

bool Output_Pair(Alignment A1,Alignment A1_P,Alignment B1,Alignment B1_P,int Read_Length,int Read_Length2)
{
	int CUTOFF=int(SW_SIMILARITY_FOR_RESCUE*Read_Length*match/100);
	if(LENGTH_CUTOFF)
	{
		CUTOFF=LENGTH_CUTOFF*match;
	}
	if(labs(A1.Loc-B1.Loc)<Read_Length)// && labs(A1_P.Loc-B1_P.Loc)>Read_Length)
	{
		if(B1_P.SW_Score> CUTOFF)
		{
			return true;
		}
		if(A1_P.SW_Score< CUTOFF)
		{
			return true;
		}
	}
	if(labs(A1_P.Loc-B1_P.Loc)<Read_Length2)// && labs(A1_P.Loc-B1_P.Loc)>Read_Length)
	{
		if(B1.SW_Score> CUTOFF)
		{
			return true;
		}
		if(A1.SW_Score< CUTOFF)
		{
			return true;
		}
	}
	else
	{
		if(B1.SW_Score> CUTOFF && B1_P.SW_Score>CUTOFF)
		{
			return true;
		}
	}
	return false;
}

//TODO: Put this into mismatch finder..
void Fix_Offset(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off)
{
	while(!A.empty())
	{
		Alignment Aln=A.top();A.pop();
		if(Aln.Sign == '+')
			Aln.Loc-=Offset;
		else
		{
			if(Aln.Sign != '-') continue;
			Aln.Loc+=Neg_Off;
		}
		if(Aln.Loc>=0)
		{
			T.push(Aln);
		}
	}
}
void Fix_Offset_Tail(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off){
        while(!A.empty())
        {    
                Alignment Aln=A.top();A.pop();
                if(Aln.Sign == '+')
                        Aln.Loc-=Offset;

                if(Aln.Loc>=0)
                {    
                        T.push(Aln);
                }            
        }            
}

void Fix_Offset_Head(std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & A,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & T,int Offset,int Neg_Off){
        while(!A.empty())
        {    
                Alignment Aln=A.top();A.pop();
                if(Aln.Sign == '-')
                        Aln.Loc-=Neg_Off;

		if(Aln.Loc>=0)
                {    
                        T.push(Aln);
                }            
        }            
}
void Print_Pair(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,Hit_Info &  HitsH,Hit_Info &  HitsT,char Sign1,char Sign2,READ & RawR,READ & RawM,unsigned char* Original_Text,char source1,char source2,FILE* Single_File,Final_Hit & H,Final_Hit & T,READ & R1, READ & R2,MEMX & MF,MEMX & MC,BWT* fwfmi,BWT* revfmi)
{
	bool Bad_Hit=false;
	int Insert_Size;
	unsigned Proper_Pair=0;
	int HP_Clip,HS_Clip,TP_Clip,TS_Clip;
	if(DASH_DEL)
	{
		if(R1.Description[strlen(R1.Description)-2]=='/')
		{
			R1.Description[strlen(R1.Description)-2]=0;
		}
		if(R2.Description[strlen(R2.Description)-2]=='/')
		{
			R2.Description[strlen(R2.Description)-2]=0;
		}
	}

	if (H.Loc != INT_MAX && T.Loc != INT_MAX && false)
	{
		Ann_Info Ann1,Ann2;
		H.Skip=Get_Skip(H.CIG);
		Location_To_Genome(H.Loc,Ann1);
		H.Flag|=R1.Tag_Number;H.Flag|=PE_SEQ;
		if (H.Loc+H.Skip > Ann1.Size) 
		{
			Bad_Hit=true;
			H.Flag |= 4;
		}
		Rescue_Clip(H.CIG,HP_Clip,HS_Clip);//Calculate clip size..
		
		T.Skip=Get_Skip(T.CIG);
		Location_To_Genome(T.Loc,Ann2);
		T.Flag|=R2.Tag_Number;T.Flag|=PE_SEQ;
		if (T.Loc+T.Skip > Ann2.Size) 
		{
			Bad_Hit=true;
			T.Flag |= 4;
		}
		Rescue_Clip(T.CIG,TP_Clip,TS_Clip);

		bool Same_Chrome=(!strcmp(Ann1.Name,Ann2.Name))? true:false;

		Insert_Size=0;
		if(H.Flag & 16)
		{
			T.Flag|=MATE_MINUS;
			if(!(T.Flag & 16))//5335 orientation?
			{
				if(Same_Chrome)
				{
					Insert_Size= -(H.Loc-T.Loc+T.Skip);
					if(Check_Proper_Pair('-','+',H.Loc,T.Loc,0/*R1.Real_Len*/))
					{
						Proper_Pair=PROPER_PAIR;
					}
				}
			}
			else
			{
				H.Flag|=MATE_MINUS;
			}
		}
		else
		{
			if(T.Flag & 16)//5335 orientation?
			{
				H.Flag|=MATE_MINUS;
				if(Same_Chrome)
				{
					Insert_Size=-H.Loc+T.Loc+T.Skip;
					if(Check_Proper_Pair('+','-',H.Loc,T.Loc,0/*R1.Real_Len*/))
					{
						Proper_Pair=PROPER_PAIR;
					}
				}
			}
		}
		if(ESTIMATE)
		{
			if(Insert_Size)
			{
				pthread_mutex_lock(&Lock_Estimate);
				Estimate.push_back((Insert_Size>0)?Insert_Size:-Insert_Size);	
				pthread_mutex_unlock(&Lock_Estimate);
				return;
			}
			else
			{
				return;
			}
		}
		H.Flag|=Proper_Pair;T.Flag|=Proper_Pair;

		//int ed1=Get_ED(H.CIG),ed2=Get_ED(T.CIG);
		//if( H.Mismatch <=BP.MAX_MISMATCHES && T.Mismatch <=BP.MAX_MISMATCHES) //((R1.Real_Len*0.02+1 >= ed1) || (R2.Real_Len*0.02+1 >= ed2) ) && //MISMATCH
		{
			Alignment_Pair Pai; 
			
		    if(H.hitType == '1' || H.hitType == '4') {
			Pai.chrom1=Ann1.Name;Pai.Loc1=H.Loc;Pai.Flag1=H.Flag;Pai.source1=H.hitType;Pai.Mismatch1=H.Mismatch;Pai.ReadLen1=R1.Real_Len;Pai.Cigar1=H.CIG;
			strcpy(Pai.Description1,R1.Description);strcpy(Pai.Tag_Copy1,RawR.Tag_Copy);strcpy(Pai.Quality1,RawR.Quality);
			
			Pai.chrom2=Ann2.Name;Pai.Loc2=T.Loc;Pai.Flag2=T.Flag;Pai.source2=T.hitType;Pai.Mismatch2=T.Mismatch;Pai.ReadLen2=R2.Real_Len;Pai.Cigar2=T.CIG;
			strcpy(Pai.Description2,R2.Description);strcpy(Pai.Tag_Copy2,RawM.Tag_Copy);strcpy(Pai.Quality2,RawM.Quality);
			Pai.Q1=H.Quality_Score;Pai.Q2=T.Quality_Score;//printf("\n=== %d %d || %d %d || %d %d\n",H.Score,T.Score,H.QScore,T.QScore,H.Quality_Score,T.Quality_Score);
			Pai.Insert_Size=Insert_Size;
		    }else {
			Pai.chrom2=Ann1.Name;Pai.Loc2=H.Loc;Pai.Flag2=H.Flag;Pai.source2=H.hitType;Pai.Mismatch2=H.Mismatch;Pai.ReadLen2=R1.Real_Len;Pai.Cigar2=H.CIG;
                        strcpy(Pai.Description2,R1.Description);strcpy(Pai.Tag_Copy2,RawR.Tag_Copy);strcpy(Pai.Quality2,RawR.Quality);

                        Pai.chrom1=Ann2.Name;Pai.Loc1=T.Loc;Pai.Flag1=T.Flag;Pai.source1=T.hitType;Pai.Mismatch1=T.Mismatch;Pai.ReadLen1=R2.Real_Len;Pai.Cigar1=T.CIG;
                        strcpy(Pai.Description1,R2.Description);strcpy(Pai.Tag_Copy1,RawM.Tag_Copy);strcpy(Pai.Quality1,RawM.Quality);
                        Pai.Q2=H.Quality_Score;Pai.Q1=T.Quality_Score;
			Pai.Insert_Size=-Insert_Size;
		    }
			Pai.QualityScore=H.Quality_Score+T.Quality_Score;Pai.Same_Chrom=Same_Chrome;
			Pai.Mismatch=Pai.Mismatch1+Pai.Mismatch2;
			Pai.paired=0;//false;
			Pai.Score=H.Score+T.Score;Pai.Indel=H.Indel+T.Indel;
			Pai.Clip=H.Clip+T.Clip;
			Alignments_Reslut.push(Pai);


		}

	}
	else 
	{
		assert(!ESTIMATE);
		if (H.Loc == INT_MAX)
		{
			assert(T.Loc!=INT_MAX);

			Ann_Info Ann2;
			Location_To_Genome(T.Loc,Ann2);
			//T.Flag|=R2.Tag_Number;
//T.Flag|=PE_SEQ;
			T.Skip=Get_Skip(T.CIG);
			if (T.Loc+T.Skip > Ann2.Size) 
			{
				Bad_Hit=true;
				T.Flag |= 4;
			}
//			T.Flag|=MATE_UNMAPPED;
			H.Flag=4;
//H.Flag|=PE_SEQ;
//H.Flag|=R1.Tag_Number;
//			if(T.Flag & 16) H.Flag|=MATE_MINUS;
			R1.Tag_Copy[R1.Real_Len]=0;
			R1.Quality[R1.Real_Len]=0;

			if(Sign2!='-' && Sign2!='+' )
			{
				if(T.Flag & 16) Sign2='-';else Sign2='+';
			}
			int ed = Get_ED(T.CIG);
			//if( (R2.Real_Len*0.02+1 >= ed) && T.Mismatch <=BP.MAX_MISMATCHES)//T.Quality_Score>20 && !Print_hits &&
//printf("\nStore %ld %d\n", T.Loc, T.Quality_Score);
			if(T.Quality_Score<=60 && T.Quality_Score>0)
			{
				Alignment_Pair Pai; 
				Pai.chrom1=Ann2.Name;Pai.Loc1=T.Loc;Pai.Flag1=T.Flag;Pai.source1=T.hitType;Pai.Mismatch1=T.Mismatch;Pai.ReadLen1=R2.Real_Len;Pai.Cigar1=T.CIG;
				strcpy(Pai.Description1,R2.Description);strcpy(Pai.Tag_Copy1,RawM.Raw_Tag_Copy);
				strcpy(Pai.Quality1,R2.Quality);
				Pai.QualityScore=T.Quality_Score;Pai.paired=0;//false;
				Pai.Score=T.Score;Pai.Indel=T.Indel;Pai.Clip=T.Clip;
				Pai.Mismatch=T.Mismatch;
				Pai.swscore1 = T.swscore;
				Alignments_Reslut.push(Pai);

			}
		}
		if(H.Loc != INT_MAX)
		{
			assert(H.Loc!=INT_MAX);

			Ann_Info Ann1;
			Location_To_Genome(H.Loc,Ann1);//H.Chr=Ann1.Name;
			//H.Flag|=R1.Tag_Number;
//H.Flag|=PE_SEQ;
			H.Skip=Get_Skip(H.CIG);
			if (H.Loc+H.Skip > Ann1.Size) 
			{
				Bad_Hit=true;
				H.Flag |= 4;
			}
//			H.Flag|=MATE_UNMAPPED;
			T.Flag=4;
//T.Flag|=PE_SEQ;T.Flag|=R2.Tag_Number;
//			if(H.Flag & 16) T.Flag|=MATE_MINUS;
			R2.Tag_Copy[R2.Real_Len]=0;
			R2.Quality[R2.Real_Len]=0;

			if(Sign1!='-' && Sign1!='+' )
			{
				if(H.Flag & 16) Sign1='-';else Sign1='+';
			}
			int ed = Get_ED(H.CIG);
			if(H.Quality_Score<=60 && H.Quality_Score>0)
			{
				Alignment_Pair Pai; 
				Pai.chrom1=Ann1.Name;Pai.Loc1=H.Loc;Pai.Flag1=H.Flag;Pai.source1=H.hitType;Pai.Mismatch1=H.Mismatch;Pai.ReadLen1=R1.Real_Len;Pai.Cigar1=H.CIG;
				strcpy(Pai.Description1,R1.Description);strcpy(Pai.Tag_Copy1,RawR.Raw_Tag_Copy);strcpy(Pai.Quality1,R1.Quality);
				Pai.QualityScore=H.Quality_Score;Pai.paired=0;//false;
				Pai.Score=H.Score;Pai.Indel=H.Indel;Pai.Mismatch=H.Mismatch;
				Pai.Clip=H.Clip;
				Pai.swscore1 = H.swscore;
				Alignments_Reslut.push(Pai);
			}
		}
	}
}

int Get_Skip(std::string & S)
{
	int Skip=0;
	if(S.empty() ) return Skip;
	//assert(!S.empty());
	std::size_t Loc = S.find_first_of("MDIS");
	std::size_t Last_Loc=0;
	const char *C_String=S.c_str();
	assert(C_String);
	while (Loc!=std::string::npos)
	{
		if((C_String[Loc] == 'M') || (C_String[Loc] == 'I'))
		{
			Skip+=atoi(C_String+Last_Loc);
		}
		Last_Loc=Loc+1;
		Loc=S.find_first_of("MDIS",Loc+1);
	}
	return Skip;
}
int Get_ED(std::string & S)//edit distance
{
	int ed=0;
	int i;
	for(i=0;i<S.size();i++)//momomo
	{
		if(S[i]=='I' || S[i]=='D') ed++;
	}
	return ed;
}

bool Check_Proper_Pair(int S1,int S2,unsigned Loc1, unsigned Loc2,int Extra_Bit)
{
	if(S1==S2)
		return false;
	if(S1=='+')
	{
		assert(S2=='-');
		if(Loc1<=Loc2)
		{
			//if(Loc2-Loc1<=INSERTSIZE+2*STD+Extra_Bit)
			if(Loc2-Loc1<=INSERTSIZE+3*REAL_STD+Extra_Bit || Loc2-Loc1<=550)
				return true;
		}
	}
	else
	{
		assert(S2=='+');
		if(Loc1>=Loc2)
		{
			//if(Loc1-Loc2<=INSERTSIZE+2*STD+Extra_Bit)
			if(Loc1-Loc2<=INSERTSIZE+3*REAL_STD+Extra_Bit || Loc1-Loc2<= 550 )
				return true;
		}
	}
	return false;
}


void Estimate_Insert(int & INSERTSIZE,int & STD)
{
	int Size=Estimate.size();

	if(Size<15*(READS_TO_ESTIMATE/100))
	{
		fprintf(stderr,"Variability of insert size is high: Assuming wild distribution..\n");
		INSERTSIZE=800;//1000;
		REAL_STD=STD=100;//250;
		SW_STRING_BUFFER=2400;
		return;

	}

	std::sort (Estimate.begin(),Estimate.end());
	int Quarter=Size/4;
	int Q1=Quarter,Q2=2*Quarter,Q3=3*Quarter;
	int Upper=Estimate[Q3]+100;
	int Lower=Estimate[Q1]-100;
	int j=0;
	unsigned Total=0;
	for(int i=0;i<Size;i++)
	{
		//printf("%d\n",Estimate[i]);
		if(Estimate[i]>Upper)
		{
			break;
		}
		else
		{
			if(Estimate[i]>=Lower)
			{
				Estimate[j++]=Estimate[i];
				Total+=Estimate[i];
			}
		}
	}
	Estimate.resize(j);

	float Mean=Total/j;
	float N=0;
	for (int i=0;i<j;i++)
	{
		N=N+pow((Estimate[i]-Mean),2);
	}

	if(j<15*(READS_TO_ESTIMATE/100))
	{
		fprintf(stderr,"Variability of insert size is high: Assuming wild distribution..\n");
		INSERTSIZE=300;//1000
		STD=60;//250;
		SW_STRING_BUFFER=2400;
		return;
	}

	int Standard_Deviation = int(sqrt (N/j));
	REAL_STD=Standard_Deviation;
	if(Standard_Deviation>STD)
		STD=Standard_Deviation;
	INSERTSIZE=int(Mean);

	fprintf (stderr,"\rMean Insert Size detected: %u\nStandard deviation: %u/%d\n",INSERTSIZE,STD,REAL_STD);
}

void Rescue_Clip(std::string S,int & P_Clip,int & S_Clip)
{
	if(S.empty() ) return;
	assert(!S.empty());
	std::size_t Loc = S.find_first_of("MDIS");
	std::size_t Last_Loc=0;
	const char *C_String=S.c_str();
	assert(C_String);

	S_Clip=0;P_Clip=0;

	while (Loc!=std::string::npos)
	{
		if(C_String[Loc] == 'S')
		{
			if(P_Clip) 
			{
				S_Clip=atoi(C_String+Last_Loc);
			}
			else
			{
				P_Clip=atoi(C_String+Last_Loc);
			}

		}
		else
		{
			if(!P_Clip)
			{
				P_Clip=INT_MAX;
			}
		}
		Last_Loc=Loc+1;
		Loc=S.find_first_of("MDIS",Loc+1);
	}
	if(P_Clip==INT_MAX)
	{
		P_Clip=0;
	}
}

void Print_Aux_Hit(std::priority_queue <Alignment_Pair,std::vector <Alignment_Pair>,Comp_Align_Pair> & Alignments_Reslut,char Sign,READ & RawR,char source,FILE* Single_File,Final_Hit & H,Final_Hit & Aux,READ & R,int Clip1,int Clip2,MEMX & MF)
{
	Ann_Info Ann1;
	if(Aux.Score ==INT_MAX)
		return; 
	Aux.Skip=Get_Skip(Aux.CIG);
	int Clip=Clip1+Clip2;
	
	int AP_Clip,AS_Clip;
	Rescue_Clip(Aux.CIG,AP_Clip,AS_Clip);
	Location_To_Genome(Aux.Loc,Ann1);//H.Chr=Ann1.Name;
	Aux.Flag|=R.Tag_Number;Aux.Flag|=PE_SEQ;Aux.Flag|=AUX;

	if (Aux.Loc+Aux.Skip > Ann1.Size) 
	{
		return;
	}

	int C1=R.Real_Len-Clip;
	int C2=R.Real_Len-AP_Clip-AS_Clip;
	assert(C1>0 && C2>0);
	//int Covered_Percent=100*(C1+C2)/R.Real_Len;
	int Real_Cover=(R.Real_Len-(std::max(Clip1,AP_Clip))-(std::max(Clip2,AS_Clip)));
	int Covered_Percent=100*(Real_Cover)/R.Real_Len;

	if(Clip1 && Clip2)//Both ends bad..
	{
		if((R.Real_Len-AS_Clip)-Clip1>8)//Rescue and aux overlap too much..
		{
			H.Quality_Score=INT_MAX;
		}
		if((R.Real_Len-AP_Clip)-Clip2>8)//Rescue and aux overlap too much..
		{
			H.Quality_Score=INT_MAX;
		}
		if(C1 <10) //very short recovery by SW..
		{
			H.Quality_Score=INT_MAX;
		}
		//if(std::min(Clip2,Clip1)>=5)
		//	H.Quality_Score=INT_MAX;

		if(Covered_Percent >=80)
		{
			if(H.Quality_Score==0 && Aux.Quality_Score==0)
				H.Quality_Score=20;
			if(Aux.Quality_Score<H.Quality_Score)
				Aux.Quality_Score=H.Quality_Score;
			else
				H.Quality_Score=Aux.Quality_Score;
		}
		else
		{
			if(C2<20)//R.Real_Len/2)
			{
				Aux.Quality_Score=0;
			}
		}

		if(Clip1>=CLIP_SAVE_LENGTH && Clip2>=CLIP_SAVE_LENGTH)//Both sides are searched for aux hit..
		{
			H.Quality_Score=0;
			return;
		}
			
	}
	if(H.Quality_Score!=INT_MAX)	
	{
		if(Covered_Percent >=80)
		{
			if(H.Quality_Score==0 && Aux.Quality_Score==0)
				H.Quality_Score=20;
			if(Aux.Quality_Score<H.Quality_Score)
				Aux.Quality_Score=H.Quality_Score;
			else
				H.Quality_Score=Aux.Quality_Score;
		}
		else
		{
			if(C2<20)//R.Real_Len/2)
			{
				Aux.Quality_Score=0;
			}
		}
	}
	else
		H.Quality_Score=0;
	if(Aux.Flag & 16) Sign='-';else Sign='+';
/*	if( (source=='3' || source=='4'))
	{ 
		if(Sign=='+')
			Sign='-';
		else Sign=='+'; 
	}
*/	
	//int ed = Get_ED(Aux.CIG);
	//if((R.Real_Len*0.02+1 >= ed) && Aux.Mismatch <=BP.MAX_MISMATCHES) //Aux.Quality_Score>20 && 
	if(Aux.Quality_Score<=60 && Aux.Quality_Score>0)
	{
		Alignment_Pair Pai; 
		Pai.chrom1=Ann1.Name;Pai.Loc1=Aux.Loc;Pai.Flag1=Aux.Flag;Pai.source1=source;Pai.Mismatch1=Aux.Mismatch;Pai.ReadLen1=R.Real_Len;Pai.Cigar1=Aux.CIG;
		strcpy(Pai.Description1,R.Description);strcpy(Pai.Tag_Copy1,RawR.Raw_Tag_Copy);strcpy(Pai.Quality1,RawR.Quality);	
		Pai.QualityScore=Aux.Quality_Score;Pai.paired=0;//false;
		Pai.Indel=Aux.Indel;Pai.Clip=Aux.Clip;
		Pai.Mismatch=Aux.Mismatch;
		Pai.swscore1 = Aux.swscore;
		Alignments_Reslut.push(Pai);
	}
	
}

/* Given a base name, will find file names of related FM Index */
void Build_Names(const char* Genome_Name,FMFILES & F,BATPARAMETERS & BP)//build FM index file names...
{
	int Last_Dash=0;
	char* tempoptarg= (char*)Genome_Name;
	char* Name=(char*)Genome_Name;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

	char* Command_Line_Buffer;

			Command_Line_Buffer=(char*)malloc(6000);
				F.REVBWTINDEX = (char*)Command_Line_Buffer;
				//printf("%s\n",F.REVBWTINDEX);
				if(Last_Dash) Last_Dash=Genome_Name-tempoptarg+1; else Genome_Name--;
				strncpy(F.REVBWTINDEX,tempoptarg,Last_Dash);
				F.REVBWTINDEX[Last_Dash+0]='r';F.REVBWTINDEX[Last_Dash+1]='e';F.REVBWTINDEX[Last_Dash+2]='v';
				strcpy(F.REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(F.REVBWTINDEX+Last_Dash+3,".bwt"); 
//printf("%s\n",F.REVBWTINDEX);
				F.BWTFILE=F.REVBWTINDEX+500;
				strncpy(F.BWTFILE,tempoptarg,Last_Dash);
				strcpy(F.BWTFILE+Last_Dash,Genome_Name+1);
				strcat(F.BWTFILE+Last_Dash,".bwt"); 


				F.REVOCCFILE = F.BWTFILE+500;
				strncpy(F.REVOCCFILE,tempoptarg,Last_Dash);
				F.REVOCCFILE[Last_Dash+0]='r';F.REVOCCFILE[Last_Dash+1]='e';F.REVOCCFILE[Last_Dash+2]='v';
				strcpy(F.REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(F.REVOCCFILE+Last_Dash+3,".fmv"); 


				F.OCCFILE=F.REVOCCFILE+500;			
				strncpy(F.OCCFILE,tempoptarg,Last_Dash);
				strcpy(F.OCCFILE+Last_Dash,Genome_Name+1);
				strcat(F.OCCFILE+Last_Dash,".fmv"); 

				F.SAFILE=F.OCCFILE+500;			
				strncpy(F.SAFILE,tempoptarg,Last_Dash);
				strcpy(F.SAFILE+Last_Dash,Genome_Name+1);
				strcat(F.SAFILE+Last_Dash,".sa");

				F.REVSAFILE = F.SAFILE+500;
				strncpy(F.REVSAFILE,tempoptarg,Last_Dash);
				F.REVSAFILE[Last_Dash+0]='r';F.REVSAFILE[Last_Dash+1]='e';F.REVSAFILE[Last_Dash+2]='v';
				strcpy(F.REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(F.REVSAFILE+Last_Dash+3,".sa"); 

				F.BINFILE=F.REVSAFILE+500;			
				strncpy(F.BINFILE,tempoptarg,Last_Dash);
				strcpy(F.BINFILE+Last_Dash,Genome_Name+1);
				strcat(F.BINFILE+Last_Dash,".pac");

				if(!BP.USELOCATION)
				{

					F.LOCATIONFILE=F.BINFILE+500;			
					strncpy(F.LOCATIONFILE,tempoptarg,Last_Dash);
					strcpy(F.LOCATIONFILE+Last_Dash,Genome_Name+1);
					strcat(F.LOCATIONFILE+Last_Dash,".ann.location");
				}

                                F.BLKFILE = F.LOCATIONFILE+500;
                                strncpy(F.BLKFILE,tempoptarg,Last_Dash);
                                strcpy(F.BLKFILE+Last_Dash,Genome_Name+1);
                                strcat(F.BLKFILE+Last_Dash,".blk.");

                                F.INDFILE = F.BLKFILE+500;
                                strncpy(F.INDFILE,tempoptarg,Last_Dash);
                                strcpy(F.INDFILE+Last_Dash,Genome_Name+1);
                                strcat(F.INDFILE+Last_Dash,".ind.");

                                F.RANGEFILE = F.INDFILE+500;
                                strncpy(F.RANGEFILE,tempoptarg,Last_Dash);
                                strcpy(F.RANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(F.RANGEFILE+Last_Dash,".range");

                                F.NLOCATIONFILE = F.RANGEFILE+500;
                                strncpy(F.NLOCATIONFILE,tempoptarg,Last_Dash);
                                strcpy(F.NLOCATIONFILE+Last_Dash,Genome_Name+1);
                                strcat(F.NLOCATIONFILE+Last_Dash,".N.location");
}
// G-->A
inline void ReplaceGtoA(READ & R) //char* Read
{
	int i;
        for (i=0;i<R.Real_Len;i++)
        {
                if (R.Tag_Copy[i] == 'G' || R.Tag_Copy[i] == 'g') R.Tag_Copy[i]='A';
        }
}
//C-->T
inline void ReplaceCtoT(READ & R)
{
	int i;
        for (i=0;i<R.Real_Len;i++)
        {
                if (R.Tag_Copy[i]  == 'C' || R.Tag_Copy[i] == 'c') R.Tag_Copy[i]='T';
        }
}
inline unsigned char Hash(char* S)
{
	unsigned char C=0;
	for (int i=2;S[i];C+=i*S[i++]);
	return C;
}
void MEM_STRUC(MEMX & MF, MEMX & MC,MEMX & MF2,MEMX & MC2,MEMX & MCLH,MEMX & MFLH, MEMX & MFLT,MEMX & MCLT,MEMX & MFH,MEMX & MCH,MEMX & MFT,MEMX & MCT,MEMLOOK & MLook,LEN & L,MEMX & mF,MEMX & mC)
{
	MF.L=MC.L=MFLH.L=MCLH.L=MFLT.L=MCLT.L=MFH.L=MCH.L=MFT.L=MCT.L=L;
/*	int LOOKUPSIZE=3;
	MF.Lookupsize=LOOKUPSIZE;MC.Lookupsize=LOOKUPSIZE;MFLH.Lookupsize=LOOKUPSIZE;MCLH.Lookupsize=LOOKUPSIZE;
	MFLT.Lookupsize=LOOKUPSIZE;MCLT.Lookupsize=LOOKUPSIZE;
	MFH.Lookupsize=LOOKUPSIZE;MCH.Lookupsize=LOOKUPSIZE;MFT.Lookupsize=LOOKUPSIZE;MCT.Lookupsize=LOOKUPSIZE;
*/	// initialise memory structures 
	Copy_MEM(MLook,MF,MC,MAX_MISMATCHES);
	Copy_MEM(MLook,MF2,MC2,MAX_MISMATCHES);
	Copy_MEM(MLook,MFLH,MCLH,MAX_MISMATCHES);
	Copy_MEM(MLook,MFLT,MCLT,MAX_MISMATCHES);
	Copy_MEM(MLook,MFH,MCH,MAX_MISMATCHES);
	Copy_MEM(MLook,MFT,MCT,MAX_MISMATCHES);
	Copy_MEM(MLook,mF,mC,MAX_MISMATCHES);
}
