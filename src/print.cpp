#include "print.h"
#include <cstdio>
#include "stdlib.h"
#include "assert.h"
#include <math.h>
#include <algorithm>

const int SOFTCLIPLENGTH=10;
extern bool DASH_DEL;
extern int INDELGAP;
extern bool STACK_LOWQ;
int Calc_MapQ(Hit_Info & H,Alignment & A,int Clip_Count);
extern int Top_Penalty;
extern int JUMP;
extern int SCOREGAP;
extern int SW_THRESHOLD;
extern std::string RGID;
extern bool Hard_Penalty;
extern bool REALN;
extern int match_SC;

int Find_Cigar(char* Current_Tag_raw,READ & RawR,unsigned char* Original_Text,char* Cigar,Hit_Info & H,char* Current_Tag,int StringLength,READ & R,int & Clip_H,int & Clip_T)
{
	s_align* Aln;
	Ann_Info Ann;
	char Org_String[StringLength+50];
	char Org_String_Ori[StringLength+50];
	Cigar_Info Cig_Info;

	int Jump=0;if(H.Sign=='-') Jump= 0+JUMP;
	s_profile* p = ssw_init((int8_t*)Current_Tag, StringLength, mata_SC, n, 1);
	if((long)(H.Loc+Jump)<=0) H.Loc = -Jump+1;
	Get_Bases(Original_Text,H.Loc+Jump,StringLength+INDELGAP,Org_String);
	Get_Bases(Original_Text_Ori,H.Loc+Jump,StringLength+INDELGAP,Org_String_Ori);
	Aln=mengyao_ssw_core(Org_String,StringLength, Current_Tag,StringLength+INDELGAP,0,0/*DP*/, p);
	if(Aln->score1 >= ACC_SCORE)
	{
		H.SW_Score=Aln->score1;H.SW_Sub_Opt_Score=Aln->score2;
		Clip_H=Aln->read_begin1;Clip_T=0;
		if(Aln->read_end1!=StringLength-1) Clip_T=StringLength-1-Aln->read_end1;
		ssw_cigar_processQ(H.Sign,Current_Tag_raw,Org_String_Ori,Aln,Cig_Info,Org_String,Aln->ref_begin1,Current_Tag,Aln->read_begin1,StringLength,R.Quality,NULL,Clip_H,Clip_T,Hard_Penalty);
		ssw_cigar_print(H.Sign,Aln->read_begin1,Current_Tag_raw+Aln->read_begin1,Org_String_Ori+Aln->ref_begin1,Aln,Cigar,Cig_Info,Org_String+Aln->ref_begin1,Current_Tag+Aln->read_begin1,Clip_H,Clip_T,StringLength);
		H.Score= -Cig_Info.Score;
		H.QScore=Cig_Info.QScore;
		H.BQScore=Cig_Info.BQScore;
		H.Mismatch=Cig_Info.Mis;
		H.Indel=Cig_Info.Indel_Count;
		H.swscore=Cig_Info.swscore;
		//A.Sign=H.Sign;

		H.Loc=H.Loc+Jump+Aln->ref_begin1;
		align_destroy(Aln);
		init_destroy(p); 
		return (Cig_Info.Length+Clip_T);
	}
	align_destroy(Aln);
	init_destroy(p); 
	return 0;
}

bool Quality_Passed(Cigar_Info & C,int StringLength)
{
	if(C.Mis>5) return false;
	if(C.Length<StringLength)
	{
		if(StringLength-C.Length+C.Mis >5) return false;
		//if(StringLength-C.Length>10) return false;

	}
	return true;
}

void Reverse_Tag(char* Dest,const READ & R,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Char_To_CharC[R.Tag_Copy[i]];Dest++;
	}
	*Dest=0;
}

void Read2Bin(char* Dest,char* Source,int StringLength)
{
	for (int i=0;i<StringLength;i++)
	{
		*Dest=Char_To_Code[Source[i]];Dest++;
	}
	*Dest=0;
}

void Read2RevCBin(char* Dest,char* Source,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Char_To_CodeC[Source[i]];Dest++;
	}
	*Dest=0;
}

void Cigar_Check_And_Print(READ & RawR,unsigned char* Original_Text,Hit_Info &  H,BATREAD & Read,int StringLength,Final_Hit & Printed_Hit,READ & R,bool Dont_Check_Quality,int Quality_Score,Alignment A,int Clip_H,int Clip_T,char* CIG)
{
	if(H.Org_Loc!=UINT_MAX && H.Org_Loc>0)
		H.Loc=H.Org_Loc;
	Print_Sam(RawR,Original_Text,Printed_Hit,R,H,StringLength,Quality_Score,A,Clip_H,Clip_T,CIG);
}

void Print_Sam(READ & RawR,unsigned char* Original_Text,Final_Hit & Printed_Hit,READ & R,Hit_Info & H,int StringLength,int Quality_Score,Alignment A,const int TClip_H,const int TClip_T,char* TCIG)
{
	bool Skip_Penalty=true;
	int Flag=0;
	int Skip=0;//StringLength;
	char Printed_CIG[MAX_SIGLEN];
	char* Qual=R.Quality,Rev_Qual[MAXTAG];char *CIG;
	char* Tag=R.Tag_Copy,Rev_Tag[MAXTAG];
	char Rev_Tag_raw[MAXTAG];
	if(H.Loc==UINT_MAX)
	{
		printf("\nUINTMAX %s %d %s %s\n",R.Description,A.Score,A.Cigar,A.chr);
		exit(0);
	}
	assert(H.Loc!=UINT_MAX);
	if(A.Score<0) 
	{
		printf("\nSCOREMINUS %s %d %s %s\n",R.Description,A.Score,A.Cigar,A.chr);
		exit(0);
	}
	assert(A.Score >=0);
	H.BQScore=A.BQScore;
	Hit_Info HOld=H;bool forced_align=false;//Re-alignment forced..
	//int Real_Len=0;
	int Clip_H=TClip_H,Clip_T=TClip_T;
	int Clip_Realn_H=0,Clip_Realn_T=0;
	bool TCIG_blank=true;
	if(TCIG)
	{
		TCIG_blank=false;
		if(TCIG[0])
		{
			CIG=TCIG;
			if(A.Loc) {
				H.SW_Score=A.SW_Score;
				H.swscore = A.swscore;
			}
//			H.SW_Sub_Opt_Score=0;
			HOld=H;
			if(REALN)
			{
				TCIG=NULL;
				forced_align=true;
			}
		}
		else/*Should be bad Cigar*/
		{
			CIG=Printed_CIG;
			TCIG=NULL;
		}
	}
	else
	{
		assert(false);
		CIG=Printed_CIG;
	}

	int Sub_Opt_Score=0;
//printf("\n ==== %s %d %d %c %d %d %d %d %d\n",TCIG,A.Loc,A.BQScore,H.Sign,R.Real_Len,StringLength,Quality_Score, H.Mismatch, H.Loc);
	if(R.Real_Len>=StringLength)
	{
		R.Tag_Copy[R.Real_Len]=0;
		RawR.Tag_Copy[R.Real_Len]=0;
		R.Quality[R.Real_Len]=0;
		char Real_String[R.Real_Len];
		char Real_String_raw[R.Real_Len];
		if(H.Sign=='+')
		{
			if(!TCIG)
			{
				Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
				Read2Bin(Real_String_raw,RawR.Tag_Copy,R.Real_Len);
				Skip=Find_Cigar(Real_String_raw,RawR,Original_Text,CIG,H,Real_String,R.Real_Len,R,Clip_H,Clip_T);
				if(forced_align)
				{
					//H.SW_Score=HOld.SW_Score;
					H.SW_Sub_Opt_Score=HOld.SW_Sub_Opt_Score;
					if(!Hard_Penalty && (Clip_T+Clip_H>0))
					{
						Skip_Penalty=false;
					}
					else
					{
						H.BQScore=HOld.BQScore;
					}
					Clip_Realn_T=Clip_T;Clip_Realn_H=Clip_H;
					Clip_H=TClip_H;Clip_T=TClip_T;
					H.Score= HOld.Score;
					H.QScore=HOld.QScore;
					H.swscore = HOld.swscore;
				}
			}

		}
		else
		{
			Reverse_Quality(Rev_Qual,R,R.Real_Len);
			Reverse_Tag(Rev_Tag,R,R.Real_Len);
			Reverse_Tag(Rev_Tag_raw,RawR,R.Real_Len);
			for(int i=0;i<R.Real_Len;i++)
			{
				R.Tag_Copy[i]=Rev_Tag[i];
				R.Quality[i]=Rev_Qual[i];
				RawR.Tag_Copy[i]=Rev_Tag_raw[i];
			}
			if(!TCIG)
			{
				Read2Bin(Real_String,R.Tag_Copy,R.Real_Len);
				Read2Bin(Real_String_raw,RawR.Tag_Copy,R.Real_Len);
				H.Loc-=(R.Real_Len-StringLength)+INDELGAP-1;
				Skip=Find_Cigar(Real_String_raw,RawR,Original_Text,CIG,H,Real_String,R.Real_Len,R,Clip_H,Clip_T);
				A.Mismatch=H.Mismatch;
				if(forced_align)
				{
					//H.SW_Score=HOld.SW_Score;
					H.SW_Sub_Opt_Score=HOld.SW_Sub_Opt_Score;
					if(!Hard_Penalty && (Clip_T+Clip_H>0))
					{
						Skip_Penalty=true;
					}
					else
					{
						H.BQScore=HOld.BQScore;
					}
					Clip_Realn_T=Clip_T;Clip_Realn_H=Clip_H;
					Clip_H=TClip_H;Clip_T=TClip_T;
					H.Score= HOld.Score;
					H.QScore=HOld.QScore;
					H.swscore=HOld.swscore;
				}
			}
		}
		if(TCIG || !TCIG_blank)
		{
			if(A.Loc)
			{
				H.Score= A.Score;
				H.QScore=A.QScore;
				if(Skip_Penalty) // Do not penalise soft clipped..
				{
					H.BQScore=A.BQScore;
				}
				H.Mismatch=A.Mismatch;
				H.Indel=A.Indel;
				H.swscore=A.swscore;
			}
		}
		else
			H.Score= -H.Score;
	}
	else
	{
		assert(false);
		sprintf(CIG,"%dM",StringLength);
		R.Tag_Copy[StringLength]=0;
		R.Quality[StringLength]=0;
	}

	if(Quality_Score)
	{
if(H.BQScore==INT_MAX)
{
        printf("\n%s %s\n",R.Description,CIG);
        exit(0);
}
//printf("\nTTTT %d %d %d\n", Quality_Score, H.BQScore, H.Score);
		Quality_Score=Calc_MapQ(H,A,Clip_H+Clip_T);
		if(Quality_Score==1)
		{
			int SW_THRESHOLD;
			if(Clip_Realn_H+Clip_Realn_T)//Realignment clips ..
			{
				assert(R.Real_Len-Clip_Realn_T-Clip_Realn_H>0);
				SW_THRESHOLD=80*(R.Real_Len-Clip_Realn_H-Clip_Realn_T)*match_SC/100;
				A.SW_Score=H.SW_Score;
				A.swscore = H.swscore;
			}
			else
			{
				SW_THRESHOLD=80*(R.Real_Len-Clip_H-Clip_T)*match/100;
			}

			if(A.SW_Score < SW_THRESHOLD)
				Quality_Score=0;
			if(R.NCount >std::max(NCOUNT,int(5*StringLength/100)))
				Quality_Score=0;
		}
	}

	if(!CIG[0])
	{
		sprintf(CIG,"%dX",R.Real_Len);
		Quality_Score=0;
		fprintf (stdout,"\nCigar Error:%s\n",R.Description);
	}
	if(H.Sign=='-') 
	{
		Flag=16; 
		Qual=Rev_Qual;
		Tag=Rev_Tag;
	}
	else
	{
		Flag=0;
	}
//printf("\nKKKLKLKL %d %d\n", A.Sub_Opt_Score, H.SW_Sub_Opt_Score);
	if(Sub_Opt_Score!=INT_MAX)
	{
		Printed_Hit.Indel=H.Indel;
		Printed_Hit.Loc=H.Loc;
		Printed_Hit.Skip=Skip;
		Printed_Hit.Flag=Flag;
		Printed_Hit.Quality_Score=Quality_Score;
		Printed_Hit.CIG=std::string(CIG);
		Printed_Hit.Tag=std::string(Tag);
		Printed_Hit.Qual=std::string(Qual);
		Printed_Hit.Mismatch=H.Mismatch;
		Printed_Hit.Score=H.Score;
		Printed_Hit.Sub_Opt_Score=A.Sub_Opt_Score;///////
		Printed_Hit.QScore=H.QScore;
		Printed_Hit.SW_Score=H.SW_Score;
		Printed_Hit.swscore=H.swscore;
		if(Quality_Score == 0)
			Printed_Hit.SW_Sub_Opt_Score = H.SW_Score;
		else
			Printed_Hit.SW_Sub_Opt_Score=H.SW_Sub_Opt_Score;
		//
		Printed_Hit.Clip=Clip_H+Clip_T;
		Printed_Hit.hitType=H.hitType;
	}
	else
	{
		Printed_Hit.Indel=H.Indel;
		Printed_Hit.Loc=H.Loc;
		Printed_Hit.Flag=Flag;
		Printed_Hit.Skip=Skip;
		Printed_Hit.Quality_Score=Quality_Score;
		Printed_Hit.CIG=std::string(CIG);
		Printed_Hit.Tag=std::string(Tag);
		Printed_Hit.Qual=std::string(Qual);
		Printed_Hit.Mismatch=H.Mismatch;
		Printed_Hit.Score=H.Score;
		Printed_Hit.Sub_Opt_Score=INT_MAX;///////
		Printed_Hit.QScore=H.QScore;
		Printed_Hit.SW_Score=H.SW_Score;
		Printed_Hit.swscore=H.swscore;
		Printed_Hit.SW_Sub_Opt_Score=H.SW_Sub_Opt_Score;
		//
		Printed_Hit.Clip=Clip_H+Clip_T;
		Printed_Hit.hitType=H.hitType;
	}
}

/*void Print_Mishit(READ & R,FILE* Mishit_File)
{
	if(PRINT_MISHIT) fprintf(Mishit_File,"@%s\n%s\n+\n%s\n",R.Description+1,R.Tag_Copy,R.Quality);	
}

void Print_Blank_Line(FILE* Single_File,READ & R)
{
	R.Tag_Copy[R.Real_Len]=R.Quality[R.Real_Len]=0;
	fprintf(Single_File,"%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",R.Description+1,R.Tag_Copy,R.Quality);
}*/

void get_len_m(const char* Cigar, int & lenM){
        const char* cig=Cigar;
        char temp[8];
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
                        cig++;n=0;
                }else if(*cig=='M')
                {
                        temp[n]='\0';int length=atoi(temp);
                        lenM += length;
                        cig++;n=0;
                }else if(*cig=='I')
                {
                        int i;temp[n]='\0';int length=atoi(temp);
			lenM += length;
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
int Mul_mapQ(Hit_Info & H, Alignment & A) {
	int l = 0;
	get_len_m(A.Cigar, l);
        double identity = 1. - (double)(l * match - H.SW_Score) / (match + mismatch) / l;
   	int mapq = 0;
        double tmp = 0, mapQ_coef_len = 50;
	float mapQ_coef_fac=log(mapQ_coef_len);
        tmp = l < mapQ_coef_len? 1. : mapQ_coef_fac / log(l);
        tmp *= identity * identity;
        mapq = (int)(6.02 * (H.SW_Score - H.SW_Sub_Opt_Score) / match * tmp * tmp + .499);
        mapq -= (int)(4.343 * log(1/*a->sub_n*/+1) + .499);
        if (mapq > 60) mapq = 60;
        if (mapq < 0) mapq = 0;
	float frac_rep = 0;
	if(A.Mismatch > 0 && H.SW_Score - H.SW_Sub_Opt_Score <= 30) frac_rep = 1./A.Mismatch < 1. - 1./A.Mismatch?1./A.Mismatch : 1. - 1./A.Mismatch;
        mapq = (int)(mapq * (1. - frac_rep) + .499);
//fprintf(stderr, "\n%d %f %f %f\n", l, identity, mapQ_coef_len, mapQ_coef_fac);
//fprintf(stderr, "\n %d\n", mapq);
        return mapq;
}

int Cal_mapQ(int alignsize,int subsw, int swscore, const char* Cigar, int Mismatch, int readlen) {
        int l = 0;
        get_len_m(Cigar, l);
        double identity = 1. - (double)(l * match - swscore) / (match + mismatch) / l;
        int mapq = 0;
        double tmp = 0, mapQ_coef_len = 50;
        float mapQ_coef_fac=log(mapQ_coef_len);
        tmp = l < mapQ_coef_len? 1. : mapQ_coef_fac / log(l);
        tmp *= identity * identity;
        mapq = (int)(6.02 * (swscore - subsw) / match * tmp * tmp + .499);

    //    mapq -= (int)(4.343 * log(1/*a->sub_n*/+1) + .499);
    	int minusscore = alignsize > 1 ? alignsize: 1;
	mapq -= (int)(4.343 * log(minusscore+1) + .499);
        if (mapq > 60) mapq = 60;
        if (mapq < 0) mapq = 0;
        float frac_rep = 0, frac = 1. ;//l*1./readlen;
        if(Mismatch > 0 && swscore - subsw <= 30) frac_rep = frac/Mismatch < frac - frac/Mismatch?frac/Mismatch : frac - frac/Mismatch;
//	printf("\n%f %d %d %d %d\n", frac_rep, swscore, subsw, alignsize,  mapq);
        mapq = (int)(mapq * (1. - frac_rep) + .499);
        return mapq;
}

int Cal_mapQ(int subsw, int swscore, const char* Cigar, int Mismatch, int readlen) {
        int l = 0;
        get_len_m(Cigar, l); 
        double identity = 1. - (double)(l * match - swscore) / (match + mismatch) / l;
        int mapq = 0;
        double tmp = 0, mapQ_coef_len = 50; 
        float mapQ_coef_fac=log(mapQ_coef_len);
        tmp = l < mapQ_coef_len? 1. : mapQ_coef_fac / log(l);
        tmp *= identity * identity;
        mapq = (int)(6.02 * (swscore - subsw) / match * tmp * tmp + .499);

        mapq -= (int)(4.343 * log(1/*a->sub_n*/+1) + .499);
        if (mapq > 60) mapq = 60; 
        if (mapq < 0) mapq = 0;
        float frac_rep = 0, frac = 1. ;//l*1./readlen;
        if(Mismatch > 0 && swscore - subsw <= 30) frac_rep = frac/Mismatch < frac - frac/Mismatch?frac/Mismatch : frac - frac/Mismatch;
        mapq = (int)(mapq * (1. - frac_rep) + .499);
        return mapq;
}

int Calc_MapQ(Hit_Info & H,Alignment & A,int Clip_Count)
{
	int Quality_Score;
	const int QUAL_START=60,Q_GAP=10;
	assert(H.BQScore!=INT_MAX);
	//int BOPEN=gap_open,BEXT=gap_extension;// BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
	//int BOPEN=gap_open,BEXT=gap_extension;// BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
	/*if(H.Status==PAIRED_SW)
	{
		Quality_Score=5;
	}
	else*/ if(H.Status==SW_RECOVERED || H.Status==PAIRED_SW)
	{
		int Sub_Opt_Score=A.Sub_Opt_Score;
		int MapQ=H.BQScore;
/*
		if(!Hard_Penalty)
		{
			Clip_Count=std::min(Clip_Count,11);
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}
		else if(Clip_Count>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
*/
		if(A.Score<A.Sub_Opt_Score-Q_GAP/3 )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
			if(!STACK_LOWQ && Quality_Score==1)
				Quality_Score=1;
			else
				Quality_Score+=5;
		}
		else
		{
			if(STACK_LOWQ)
				Quality_Score=1;
			else
				Quality_Score=1;
		}
		//Quality_Score=0;
		assert (Quality_Score<=60 && Quality_Score>=0);
	}
	//Unique hits..
	else if(H.Status==UNIQUEHIT)
	{
//printf("\nZZZZZZ %d %d\n", Top_Penalty, H.BQScore);
		int Sub_Opt_Score=INT_MAX;
		int MapQ=0; 
		if(H.BQScore >= -10)
			MapQ = H.BQScore;
		else
			MapQ = -10*(drand48());
		int Offset=0;
/*		if(false && !Hard_Penalty)
		{
			Clip_Count=std::min(Clip_Count,11);
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}
		else if(Clip_Count>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
*/
		Quality_Score=std::max(1,QUAL_START-Offset+MapQ-std::min(Top_Penalty,QUAL_START/3));
		if(Quality_Score>1) 
		{
			Quality_Score+=5;
			if(Quality_Score>60) Quality_Score=60;
		}
		//if(!STACK_LOWQ && Quality_Score==1)
			//Quality_Score=0;
		assert (Quality_Score<=60 && Quality_Score>=0);

	}
	else if(H.Status==SHARP_UNIQUEHIT)
	{
//printf("\nGGGGGG\n");
		int Sub_Opt_Score=INT_MAX;
		int MapQ=H.BQScore;
		int Offset=0;
/*		if(false && !Hard_Penalty)
		{
			Clip_Count=std::min(Clip_Count,11);
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}
		else if(Clip_Count>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
			//Offset=5;
		}	
*/
		//if(H.Indel) MapQ+=BOPEN;
		Quality_Score=std::max(1,QUAL_START-Offset+MapQ-std::min(Top_Penalty,QUAL_START/3));
		if(Quality_Score>1) 
			Quality_Score+=5;
		//if(!STACK_LOWQ && Quality_Score==1)
		//	Quality_Score=0;
		assert (Quality_Score<=60 && Quality_Score>=0);
	}
	//Multiple hits..
	else 
	{
//printf("\nKKKK %d %d %d %d\n", A.Score, A.Sub_Opt_Score, SCOREGAP, H.BQScore);
//		int mapq = Mul_mapQ(H, A);
//		return mapq;

//printf("\nmapq %d\n", mapq);
		int Sub_Opt_Score=H.Sub_Opt_Score;
		assert(H.Status==MULTI_HIT);
		assert(H.Sub_Opt_Score!=INT_MAX);
		int MapQ=H.BQScore;
//test if should remove this
		if(!Hard_Penalty)
		{
			Clip_Count=std::min(Clip_Count,11);
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}
		else if(Clip_Count>30)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
/**/
		if(A.Score<A.Sub_Opt_Score-SCOREGAP)//Q_GAP/3 )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
			if(Quality_Score==1)
			{
				Quality_Score=std::max(1,QUAL_START/2+5+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
				if(Quality_Score==1)
				{
					Quality_Score=std::max(1,QUAL_START/2+5+5+MapQ);//-std::min(Top_Penalty,QUAL_START/2));
				}

				if(Quality_Score>5) Quality_Score=5;
				assert(Quality_Score<=5);
				assert (Quality_Score>=0);
				return Quality_Score;
			}
			else
			{
				Quality_Score+=std::max( (int) (5*(drand48()+0.2)), (int) (10*drand48()+MapQ));
				assert (Quality_Score<=60 && Quality_Score>=0);
				return Quality_Score;
			}
			/*if(!STACK_LOWQ && Quality_Score==1)
				Quality_Score=0;
			else
				Quality_Score+=5;*/
		}
		else
		{
			if(STACK_LOWQ)
				Quality_Score=1;
			else
				Quality_Score=1;
		}
		assert (Quality_Score<=60 && Quality_Score>=0);
		//Quality_Score=0;
	}
	assert (Quality_Score<=60 && Quality_Score>=0);
	return Quality_Score;
}

/*int Calc_MapQX(Hit_Info & H,Alignment & A,int Clip_Count)
{
	int Quality_Score;
	const int QUAL_START=210,Q_GAP=10;
	if(H.Status==SW_RECOVERED)
	{
		int Sub_Opt_Score=A.Sub_Opt_Score;
		int MapQ= -H.Score;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
		//MapQ= MapQ/2; 
		if(A.Score<A.Sub_Opt_Score-Q_GAP/2 )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else if(A.Score<A.Sub_Opt_Score-Q_GAP )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else
		{
			Quality_Score=0;
		}
		Quality_Score=0;
	}
	//Unique hits..
	else if(H.Status==UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ= -H.Score+H.BQScore;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ+=Score_Add;
		}	
		if(H.Indel==1) MapQ+=BOPEN;
		//MapQ= MapQ/2; 
		Quality_Score=std::max(1,QUAL_START+MapQ);

	}
	else if(H.Status==SHARP_UNIQUEHIT)
	{
		int Sub_Opt_Score=INT_MAX;
		int MapQ= -H.Score+H.BQScore;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ+=Score_Add;
		}	
		if(H.Indel==1) MapQ+=BOPEN;
		//MapQ= MapQ/2; 
		Quality_Score=std::max(1,QUAL_START+MapQ);

	}
	//Multiple hits..
	else 
	{
		int Sub_Opt_Score=H.Sub_Opt_Score;
		assert(H.Status==MULTI_HIT);
		assert(H.Sub_Opt_Score!=INT_MAX);
		int MapQ= -H.Score;
		if(Clip_Count)
		{
			int Score_Add=std::min((Clip_Count)*MN,BOPEN+BEXT*(Clip_Count-1));
			MapQ-=Score_Add;
		}	
		//MapQ= MapQ/2; 
		if(H.Score<H.Sub_Opt_Score-Q_GAP/2 )
		{
			Quality_Score=std::max(1,QUAL_START/3+MapQ);
		}
		else if(A.Score<A.Sub_Opt_Score-Q_GAP )
		{
			Quality_Score=std::max(1,QUAL_START/2+MapQ);
		}
		else
		{
			Quality_Score=0;
		}
		Quality_Score=0;
	}
	return Quality_Score;
}*/

void Print_Unmapped(FILE* Single_File,READ & R,bool Mate_Mapped,unsigned Paired,unsigned HT,int Read_Len)
{
	char* Qual=R.Quality;
	char* Tag=R.Tag_Copy;
	unsigned Flag=4;
	Flag=((Flag|Paired)|HT);
	if(!Mate_Mapped)
	{
		Flag|=8;
	}
	R.Tag_Copy[Read_Len]=R.Quality[Read_Len]=0;
	if(DASH_DEL)
	{
		if(R.Description[strlen(R.Description)-2]=='/')
		{
			R.Description[strlen(R.Description)-2]=0;
		}
	}
	fprintf(Single_File,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:%s\n",R.Description+1,Flag,R.Tag_Copy,R.Quality,RGID.c_str());
}
