#include "filters-l.h"
#include "math.h"
#include <cstdio>
#include "stdlib.h"
#include "assert.h"

int g_log_n[256];
const int POWLIMIT=300;
extern int QUALITYCONVERSIONFACTOR;
const int QGAP=100;
float POW10[POWLIMIT];

int Cmp (const void * a, const void * b)
{
	  return ( *(char*)b - *(char*)a );
}

float Pr(float Q)
{
	assert((int)Q>=0 && (int)Q<POWLIMIT-1);
	return(1-POW10[(int)Q]);
	//printf("table: %f\tlib: %f\n",POW10[(int)Q],1-pow(10,-Q/10));
	//return (1-pow(10,-Q/10));
}

float Pow10(float Q)
{
	static float Max=0;
	//if(Max<Q) {Max=Q;printf("Newmax: %f\n",Max);}
	if((int)Q<0) Q=0;
	assert((int)Q>=0);
	if((int)Q<POWLIMIT-1)
		return POW10[(int)Q];
	else
		return pow(10,-float(Q)/10);
}

void Build_Pow10()
{
	for(int Q=0;Q<POWLIMIT;Q++)
	{
		POW10[Q]=(pow(10,-float(Q)/10));
	}
}

float QSumC(SARange & S,char* Q,int Mis,int StringLength)
{
	assert(S.Start && Q && Mis>=0); 

	float Q_Total=0;
	float Prob[MAXTAG];
	int L=0;
	char TChar=Q[StringLength];
	Q[StringLength]=0;

	if (INPUT_FILE_TYPE==FA)
	{
		float Q_Value=QLIMIT;//Assume uniform low quality..
		float Penalty= -10*log10(Pr(Q_Value));
		Q_Total= Penalty*(StringLength-Mis);

		for(int i=0;i<Mis;i++)
		{
			float Q_Value=QLIMIT;//assume uniform quality 0..
			float Penalty= Q_Value/3;
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			Q_Total+=Penalty;
		}
		Q[StringLength]=TChar;
		return (Q_Total);
	} 
	else
	{
		for(char* Q1=Q;*Q1!=0 && *Q1!='\n';Q1++)
		{
			float Q_Value=*Q1-QUALITYCONVERSIONFACTOR;//make quality to integer..
			assert(Q_Value>=0);
			float Penalty= -10*log10(Pr(Q_Value));
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			//Penalty=0;
			Prob[L++]= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)
		}
		//assert(L>0);

		for(int i=0;i<Mis;i++)
		{
			float Q_Value=Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR;//make quality to integer..
			float Penalty= Q_Value/3;
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			Prob[S.Mismatch_Pos[i]]= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)
		}
		for(int i=0;i<L;i++) Q_Total+=Prob[i];
		Q[StringLength]=TChar;
		return (Q_Total);
	}

}
bool SNP_Check(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC)
{
	char Rev_Qual[MAXTAG];
	assert (false);
		/*
	memcpy(Rev_Qual,R.Quality,StringLength);
	qsort(Rev_Qual,StringLength,sizeof(char),Cmp);
	int QUAL_FOR_SNP=Rev_Qual[10]-QUALITYCONVERSIONFACTOR;

	Reverse_Quality(Rev_Qual,R,StringLength);
	if(Top_Mis<=1) return true;
	char* Q;
	SARange S;
	if(Plus_Hits)
	{
		Q=R.Quality;
		S=MF.Hit_Array[0];
	}
	else
	{
		Q=Rev_Qual;
		S=MC.Hit_Array[0];
	}

	bool Go_Deeper=false;int High_Mis_Count=0;
	//const int QUAL_FOR_SNP=20;
	for(int i=0;i<Top_Mis;i++)
	{
		if((Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR)>=QUAL_FOR_SNP) {Go_Deeper=true;High_Mis_Count++;}//make quality to integer..
	}
	if(!Go_Deeper) return true;//true;//else return false;
	if(High_Mis_Count>1) return true;//true;//else return false;

	SARange* Sub_F=MF.Hit_Array+Plus_Hits+1;
	SARange* Sub_C=MC.Hit_Array+Minus_Hits+1;

	bool All_Mismatch_Low_Quality=false;
	int Opt_In_Plus=0,Opt_In_Minus=0;
	Q=R.Quality;
	for(int i=0;Sub_F[i].Start;i++)
	{
		S=Sub_F[i];
		bool High_Quality_Base=false;
		for(int i=0;i<Top_Mis;i++)
		{
			if((Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR)>=QUAL_FOR_SNP) {High_Quality_Base=true;break;}//make quality to integer..
		}
		if(!High_Quality_Base) 
		{
			All_Mismatch_Low_Quality=true;
			Opt_In_Plus++;
			MF.Hit_Array[0]=S;
		}
	}

	Q=Rev_Qual;
	for(int i=0;Sub_C[i].Start;i++)
	{
		S=Sub_C[i];
		bool High_Quality_Base=false;
		for(int i=0;i<Top_Mis;i++)
		{
			if((Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR)>=QUAL_FOR_SNP) {High_Quality_Base=true;break;}//make quality to integer..
		}
		if(!High_Quality_Base) 
		{
			All_Mismatch_Low_Quality=true;
			Opt_In_Minus++;
			MC.Hit_Array[0]=S;
		}
	}

	if (Opt_In_Minus+Opt_In_Plus==1)
	{
		if(Opt_In_Plus) 
		{
			Plus_Hits=1;Minus_Hits=0;
			MF.Hit_Array_Ptr=2;MF.Hit_Array[1].Start=0;
			MC.Hit_Array_Ptr=1;MC.Hit_Array[0].Start=0;
		}
		else 
		{
			Plus_Hits=0;Minus_Hits=1;
			MC.Hit_Array_Ptr=2;MC.Hit_Array[1].Start=0;
			MF.Hit_Array_Ptr=1;MF.Hit_Array[0].Start=0;
		}
		return true;
	}
	else if(All_Mismatch_Low_Quality) return false;
	*/
}

bool Phred_Check(int & Plus_Hits,int & Minus_Hits,int Top_Mis,int Subopt_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,float & Top_BQScore, float & Sub_Score, float & Sub_BQScore,int & Quality_Score)
{
	SARange *Sub_F,*Sub_C,Opt_SA,Top_SA;
	char Rev_Qual[MAXTAG];

	Sub_F=MF.Hit_Array+Plus_Hits+1;
	Sub_C=MC.Hit_Array+Minus_Hits+1;
	Reverse_Quality(Rev_Qual,R,StringLength);

	float TS,Sub_Sub_Score=FLT_MAX,Prob_Sum=0;
	Sub_Score=Sub_BQScore=FLT_MAX;

	if(Plus_Hits)
	{
		Top_SA=MF.Hit_Array[0];
		Top_Score=QSumX(MF.Hit_Array[0],R.Quality,Top_Mis,StringLength);
		Top_BQScore=BQSumX(MF.Hit_Array[0],R.Quality,Top_Mis,StringLength);
		TS=QSumC(MF.Hit_Array[0],R.Quality,Top_Mis,StringLength);
	}
	else
	{
		Top_SA=MC.Hit_Array[0];
		Top_Score=QSumX(MC.Hit_Array[0],Rev_Qual,Top_Mis,StringLength);
		Top_BQScore=BQSumX(MC.Hit_Array[0],Rev_Qual,Top_Mis,StringLength);
		TS=QSumC(MC.Hit_Array[0],Rev_Qual,Top_Mis,StringLength);
	}
	
	int Opt_In_Plus=0,Opt_In_Minus=0;
	int Hits=0;
	for(int i=0;Sub_F[i].Start;i++)
	{
		assert(Sub_F[i].End-Sub_F[i].Start>=0);
		float S=QSumX(Sub_F[i],R.Quality,Subopt_Mis,StringLength),S_1=QSumC(Sub_F[i],R.Quality,Subopt_Mis,StringLength);
		float SBQ=BQSumX(Sub_F[i],R.Quality,Subopt_Mis,StringLength);
		//S_1=pow(10,-S_1/10);
		S_1=Pow10(S_1);
		Prob_Sum+=(Sub_F[i].End-Sub_F[i].Start+1)*S_1;
		Hits+=(Sub_F[i].End-Sub_F[i].Start);
		if(S<=Sub_Score) 
		{
			Sub_Sub_Score=Sub_Score;
			Sub_Score=S;
			Sub_BQScore=SBQ;
			Opt_In_Plus=1;Opt_In_Minus=0;
			Opt_SA=Sub_F[i];
		}
		else if(S<=Sub_Sub_Score)
		{
			Sub_Sub_Score=S;
		}
	}
	for(int i=0;Sub_C[i].Start;i++)
	{
		assert(Sub_C[i].End-Sub_C[i].Start>=0);
		float S=QSumX(Sub_C[i],Rev_Qual,Subopt_Mis,StringLength),S_1=QSumC(Sub_C[i],Rev_Qual,Subopt_Mis,StringLength);
		float SBQ=BQSumX(Sub_C[i],Rev_Qual,Subopt_Mis,StringLength);
		//S_1=pow(10,-S_1/10);
		S_1=Pow10(S_1);
		Prob_Sum+=(Sub_C[i].End-Sub_C[i].Start+1)*S_1;
		Hits+=(Sub_C[i].End-Sub_C[i].Start);
		if(S<=Sub_Score) 
		{
			Sub_Sub_Score=Sub_Score;
			Sub_Score=S;
			Sub_BQScore=SBQ;
			Opt_In_Plus=0;Opt_In_Minus=1;
			Opt_SA=Sub_C[i];
		}
		else if(S<=Sub_Sub_Score)
		{
			Sub_Sub_Score=S;
		}
	}
	assert(Sub_Score<=Sub_Sub_Score);
	//assert(Prob_Sum>=0 && Prob_Sum<=1);

	//float Top_Prob=pow(10,-TS/10);
	float Top_Prob=Pow10(TS);
	float Mapping_Quality,Posterior=Top_Prob/(Prob_Sum+Top_Prob);
	Mapping_Quality= std::max(0,int(-10*log10(1-Posterior)));//avoid minus infinity..
	if(LAZY_PHRED) Mapping_Quality=bwa_approx_mapQ(Hits,1);
	Quality_Score=std::min(40,int(Mapping_Quality));
	//assert(Quality_Score>=0 && Quality_Score!=30);

	if(Top_Score<Sub_Score-QGAP)// Least mismatch hit is the best..
	{
		Quality_Score=30;
		return true;
	}
	else 
	{
		return false;
	}
}

/*float QSum(SARange & S,char* Q,int Mis)
{
	assert(S.Start && Q && Mis>=0); 
	if (!Mis) return 0;

	float Q_Total=0;
	float Prob[MAXTAG];
	int L=0;

	for(char* Q1=Q;*Q1!=0 && *Q1!='\n';Q1++)
	{
		float Q_Value=*Q1-QUALITYCONVERSIONFACTOR;//make quality to integer..
		assert(Q_Value>=0);
		float Penalty=10*log10(1-pow(10,-Q_Value/10));
		Prob[L++]= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)
	}
	assert(L>0);

	for(int i=0;i<Mis;i++)
	{
		float Q_Value=Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR;//make quality to integer..
		float Penalty=10*log10((1-pow(10,-Q_Value/10))/3);
		Prob[S.Mismatch_Pos[i]]= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)
	}

	for(int i=0;i<L;i++) Q_Total+=Prob[i];
	//assert(Q_Total>0);
	return -Q_Total;
}*/

void Reverse_Quality(char* Dest,const READ & R,int StringLength)
{
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=R.Quality[i];Dest++;
	}
	*Dest=0;
}

bool Phred_Check_Multi(int & Plus_Hits,int & Minus_Hits,int Top_Mis,READ & R, int StringLength,MEMX & MF,MEMX & MC,float & Top_Score,int & Quality_Score,char & Status)
{
	char Rev_Qual[MAXTAG];
	char Opt_Sign;
	int Hits=0;
	SARange Opt_SA;

	Reverse_Quality(Rev_Qual,R,StringLength);
	Top_Score=INT_MAX;
	float Second_Score=INT_MAX,Prob_Sum=0,S_1=0;
	float TS;

	if(Plus_Hits)
	{
		for(int i=0;MF.Hit_Array[i].Start;i++)
		{
			SARange SA=MF.Hit_Array[i];
			assert(SA.End-SA.Start>=0);
			Hits+=(SA.End-SA.Start+1);
			float S=QSumX(SA,R.Quality,Top_Mis,StringLength);//-Subopt_Mis*mismatch;//+(StringLength-Subopt_Mis)*match;;	
			
			S_1=QSumC(SA,R.Quality,Top_Mis,StringLength);
			//S_1=pow(10,-S_1/10);
			//Prob_Sum+=(SA.End-SA.Start+1)*pow(10,-S_1/10);
			Prob_Sum+=(SA.End-SA.Start+1)*Pow10(S_1);
			
			if(S<=Top_Score) 
			{
				Second_Score=Top_Score;
				Top_Score=S;
				TS=S_1;
				Opt_SA=SA;Opt_Sign='+';
				assert(Top_Score<=Second_Score);
			}
			else if(S<=Second_Score) 
			{
				Second_Score=S;
				assert(Top_Score<=Second_Score);
			}
		}
	}
	if(Minus_Hits)
	{
		for(int i=0;MC.Hit_Array[i].Start;i++)
		{
			SARange SA=MC.Hit_Array[i];
			Hits+=(SA.End-SA.Start+1);
			assert(SA.End-SA.Start>=0);
			float S=QSumX(SA,Rev_Qual,Top_Mis,StringLength);

			S_1=QSumC(SA,Rev_Qual,Top_Mis,StringLength);
			//S_1=pow(10,-S_1/10);
			//Prob_Sum+=(SA.End-SA.Start+1)*pow(10,-S_1/10);
			Prob_Sum+=(SA.End-SA.Start+1)*Pow10(S_1);
			//S_1=pow(10,-S_1/10);
			//Prob_Sum+=(SA.End-SA.Start+1)*S_1;

			if(S<=Top_Score) 
			{
				Second_Score=Top_Score;
				Top_Score=S;
				TS=S_1;
				Opt_SA=SA;Opt_Sign='-';
				assert(Top_Score<=Second_Score);
			}
			else if(S<=Second_Score) 
			{
				Second_Score=S;
				assert(Top_Score<=Second_Score);
			}
		}
	}

	//float Top_Prob=pow(10,-Top_Score/10);
	//float Top_Prob=pow(10,-TS/10);
	float Top_Prob=Pow10(TS);
	float Mapping_Quality,Posterior=Top_Prob/(Prob_Sum);
	Mapping_Quality= std::max(0,int(-10*log10(1-Posterior)));//avoid minus infinity..
	if(LAZY_PHRED) Mapping_Quality=bwa_approx_mapQ(Hits,1);
	Quality_Score=std::min(40,int(Mapping_Quality));
	//assert(Quality_Score>=0 && Quality_Score!=30);


	if(Opt_SA.Start==Opt_SA.End)
	{
		if(Top_Score<Second_Score-QGAP)//use this to change accuracy.. 
		{
			if(Opt_Sign=='+') 
			{
				Plus_Hits=1;Minus_Hits=0;
				MF.Hit_Array[0]=Opt_SA;
				MF.Hit_Array_Ptr=2;MF.Hit_Array[1].Start=0;
				MC.Hit_Array_Ptr=1;MC.Hit_Array[0].Start=0;
			}
			else 
			{
				Plus_Hits=0;Minus_Hits=1;
				MC.Hit_Array[0]=Opt_SA;
				MC.Hit_Array_Ptr=2;MC.Hit_Array[1].Start=0;
				MF.Hit_Array_Ptr=1;MF.Hit_Array[0].Start=0;
			}
			return true;
		}
		else 
		{
			Status=DO_SW_SCAN;
			//Quality_Score=0;
			return true; 
		}
	}
	//Quality_Score=0;
	Status=DO_SW_SCAN;
	return true;
}

float QSumX(SARange & S,char* Q,int Mis,int StringLength)
{
	if(!(S.Start && Q && Mis>=0)) return 0;
	assert(S.Start && Q && Mis>=0); 
	//if (!Mis) return 0;

	float Q_Total=0;
	float Prob[MAXTAG];
	int L=0;
	char TChar=Q[StringLength];
	Q[StringLength]=0;

	if (INPUT_FILE_TYPE==FA)
	{
		float Q_Value=QLIMIT;//Assume uniform low quality..
		float Penalty= -10*log10(Pr(Q_Value));
		Q_Total= Penalty*(StringLength-Mis);

		for(int i=0;i<Mis;i++)
		{
			float Q_Value=QLIMIT;
			float Penalty= -10*log10((1-Pr(Q_Value))/3);
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			Q_Total+=Penalty;
		}
		Q[StringLength]=TChar;
		return Q_Total;
	} 
	else
	{
		for(char* Q1=Q;*Q1!=0 && *Q1!='\n';Q1++)
		{
			float Q_Value=*Q1-QUALITYCONVERSIONFACTOR;//make quality to integer..
			Q_Value=std::max(30.0f,Q_Value);
			assert(Q_Value>=0);// && Q_Value<=40);
			float Penalty= -10*log10(Pr(Q_Value));
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			assert(Penalty<=QLIMIT); 
			Prob[L++]= Penalty;
		}
		//assert(L>0);

		for(int i=0;i<Mis;i++)
		{
			assert(i<5 && Mis>0);
			float Q_Value=Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR;//make quality to integer..
			Q_Value=std::max(0.0f,Q_Value);
			assert(Q_Value>=0);// && Q_Value<=40);
			float Penalty= -10*log10((1-Pr(Q_Value))/3);
			Penalty=std::min(QLIMIT_FLOAT,Penalty);
			assert(Penalty<=QLIMIT); 
			Prob[S.Mismatch_Pos[i]]= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)
		}

		//for(int i=0;i<L;i++) Q_Total+=Prob[i];
		for(int i=0;i<L;i++) Q_Total+=Prob[i];
		Q[StringLength]=TChar;
		return (Q_Total);
	}
}

void bwase_initialize() 
{
	int i;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

int bwa_approx_mapQ(int p, int mm)
{
	int n;
	if (p == 0) return 37;
	n = (p >= 255)? 255 : p;
	return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}


float BQSumX(SARange & S,char* Q,int Mis,int StringLength)
{
	if(!(S.Start && Q && Mis>=0)) return 0;
	assert(S.Start && Q && Mis>=0); 
	//if (!Mis) return 0;

	float Q_Total=0;

	if (INPUT_FILE_TYPE==FA)
	{
		for(int i=0;i<Mis;i++)
		{
			float Q_Value=QLIMIT;
			Q_Total-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
		}
		return Q_Total;
	} 
	else
	{
		for(int i=0;i<Mis;i++)
		{
			float Q_Value=Q[S.Mismatch_Pos[i]]-QUALITYCONVERSIONFACTOR;//make quality to integer..
			Q_Value=std::max(0.0f,Q_Value);
			assert(Q_Value>=0);// && Q_Value<=40);
			Q_Total-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
		}
		return (Q_Total);
	}
}
