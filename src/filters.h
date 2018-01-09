#ifndef __FILTERS__
#define __FILTERS__
#include "global.h"
#include "common.h"
enum {DO_SW_SCAN,NO_SW_SCAN};

int Cmp (const void * a, const void * b);
bool SNP_Check(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC);
bool Phred_Check(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,float & Top_BQScore, float & Sub_Score, float & Sub_BQScore,int & Quality_Score);
float QSum(SARange & S,char* Q,int Mis);
void Reverse_Quality(char* Dest,const READ & R,int StringLength);
void Reverse_Quality(char* Dest,char* Quality,int StringLength);
bool Phred_Check_Multi(int & Plus_Hits,int & Minus_Hits,int Top_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,int & Quality_Score,char & Status);
float QSumX(SARange & S,char* Q,int Mis,int StringLength);
float BQSumX(SARange & S,char* Q,int Mis,int StringLength);
float QSumC(SARange & S,char* Q,int Mis,int StringLength);
void Build_Pow10();
float Pr(float Q);
float Pow10(float Q);

void bwase_initialize(); 
int bwa_approx_mapQ(int p, int mm);

#endif

