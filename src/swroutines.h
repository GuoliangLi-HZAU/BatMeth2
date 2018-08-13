#ifndef __SWROUTINES__
#define __SWROUTINES__
#include "batlib.h"
#include "rqindex.h"
#include "ssw.h"
#include "global.h"
#include "common.h"
#include <map>
#include <queue>
struct Cigar_Info
{
	int M;int I;int D;int Mis;int Length;int Score;int Indel_Count;int QScore;int BQScore;int swscore;
};
int Do_SW_Pair(char* Pattern_raw,READ & RawR,unsigned char* Original_Text,Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score,int & Total_SW_Scans);
//int Do_SW_Pair(Pair* Pairs,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int Current_Score);
int Do_SW(char* Pattern_raw,READ & RawR,RQINDEX & RQHALF,unsigned char* Original_Text,BWT* revfmi,SARange* Head_Hits,char* Pattern,int Flank_Size,int StringLength,int & Err,int Shift,char Sign,int Current_Score,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Alignments,int & Total_SW_Scans,int & filter);
void Get_Bases (unsigned char* Original_Text,unsigned Location,int StringLength,char* Org_String);
void ssw_cigar_print(char Sign,int ReadSpan,char* rawRead,char* rawRef,s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern);
void ssw_cigar_process(char Sign,int ReadSpan,char* rawRead,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,char* Pattern,int StringLength);
void ssw_cigar_processQ(char Sign,char* rawRead,char* rawRef,s_align* a,Cigar_Info & C,char* Ref,int Ref_Off,char* Pattern,int Pat_Off,int StringLength,const char* Qual,char* Cigar,int Clip_H,int Clip_T,bool Hard_Penalty=true);//,READ & R
void ssw_cigar_print(char Sign,int ReadSpan,char* rawRead,char* rawRef,s_align* a,char* Cigar,Cigar_Info & C,char* Ref,char* Pattern,int Clip_H,int Clip_T,int StringLength);
#endif
