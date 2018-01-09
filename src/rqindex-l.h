#ifndef RQINDEX_H
#define RQINDEX_H
#include "utils-l.h"

struct RQINDEX
{
	SA* SA_Index;
	int* SA_Blocks;
	char COMPRESS;
};

void Load_Range_Index(RQINDEX & R,int STRINGLENGTH,FMFILES F,unsigned & Entries);
void Load_Hits_File(BATPARAMETERS BP,FILELIST & FList,In_File & IN,char MODE);
void Search_Small_Gap(BWT* revfmi,RQINDEX & R, SARange & Head,  SARange & Tail, int d,Pair* Pairs,int & Pairs_Index,unsigned MAXCOUNT,int & HITS,unsigned Entries,unsigned Conversion_Factor);
void Get_Head_Tail(BWT* revfmi,RQINDEX & R,SARange & Head, SARange & Tail,int d,Pair* Pairs,int & Pairs_Index,unsigned MAXCOUNT,int & HITS,unsigned Entries,unsigned Conversion_Factor);
char Read_PairV(SARange* Head_Hits_Pos,SARange* Head_Hits_Neg, SARange* Tail_Hits_Pos,SARange* Tail_Hits_Neg,FILE* Data_File,In_File IN,Record_Info & Hit_Details);
char Read_Pair(BWT* fwfmi,BWT* revfmi,SARange* Head_Hits_Pos,SARange* Head_Hits_Neg, SARange* Tail_Hits_Pos,SARange* Tail_Hits_Neg,FILE* Data_File,In_File & IN,Record_Info & Hit_Details,char & MAX_HIT_FAULT);
void Load_Info( RQINDEX & R,Tag_Info & Tag, SARange & Head,unsigned Entries);
unsigned Get_Location(BWT* revfmi,RQINDEX & R,Tag_Info & Tag, unsigned Offset,unsigned Conversion_Factor);
unsigned Get_Block_Start(RQINDEX & R,unsigned SAValue,unsigned & M,unsigned H);
#endif
