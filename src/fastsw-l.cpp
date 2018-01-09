#include "fastsw-l.h"

extern int SEEDSIZE;
extern int QUALITYCONVERSIONFACTOR;
extern int QUALITYSCALEFACTOR;
extern int32_t match, mismatch, gap_open, gap_extension;
extern int gap_openP,gap_extensionP;
const int SCOREDIFF=8*match;
const bool OUTPUTALL=true;
struct Match_Info
{
	int Loc;
	int Mis;
	int Length;
};

struct Gap_Info
{
	int Start;
	int End;
	int Score;
	int Mis;
	int A,B,C;
	char Indel_Flag;
};

void Mismatch_Scan_With_ScoreX(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Ref,char* Read,int Read_Len,int Max_Mis,int Current_Mis,Alignment & A)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	int Mis=0;
	int SW_Score=0;
	float Score=0,QScore=0,BQScore=0;
	float TScore=0,TQScore=0,TBQScore=0;
	for(int j=0;j<Read_Len && Mis<=Max_Mis;j++)
	{
		if(*(Ref+j)!=*(Read+j) )//|| (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='T' && (*(rawRead+j)=='c' || *(rawRead+j)=='C' ) ) ||  (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='A' && (*(rawRead+j)=='g' ||*(rawRead+j)=='G' ) ) ) 
		{
			Mis++;
			SW_Score-=mismatch;
			Score+=Mis_FA_Score;
		}
		else
		{
			SW_Score+=match;
		}
	}
	A.Mismatch+=Mis;
	A.SW_Score+=SW_Score;
	A.Score=SW_Score;
}
void Mismatch_Scan_With_Score(char Sign,int ReadSpan,READ & RawR,char* Ref,char* Read,char* Qual,int Read_Len,int Max_Mis,int Current_Mis,Alignment & A,char* rawRef)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	int Mis=0;
	int SW_Score=0;
	float Score=0,QScore=0,BQScore=0;
	float TScore=0,TQScore=0,TBQScore=0;
	for(int j=0;j<Read_Len && Mis<=Max_Mis;j++,Qual++)
	{
		if(*(Ref+j)!=*(Read+j))// || (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='T' && (*(rawRead+j)=='c' || *(rawRead+j)=='C' ) ) ||  (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='A' && (*(rawRead+j)=='g' ||*(rawRead+j)=='G' ) )  ) 
		{
			Mis++;
			SW_Score-=mismatch;
			if(INPUT_FILE_TYPE==FA) 
			{
				Score+=Mis_FA_Score;
				QScore+=std::min(QLIMIT_FLOAT,Mis_FA_Score/3);
			}
			else
			{
				float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
				assert(Q_Value>=0);
				BQScore-= MN + floor( (MX-MN)*(std::min(Q_Value, 40.0f)/40.0f) );
				float Penalty= -10*log10((1-Pr(Q_Value))/3);
				Penalty=std::min(QLIMIT_FLOAT,Penalty);
				Score+= Penalty;//Convert to probability of base being wrong =1-10^(Q_Value/10)

				Penalty= Q_Value;///3;//prob II..
				Penalty=std::min(QLIMIT_FLOAT,Penalty/3);
				QScore+=Penalty;
			}
		}
		else
		{
			if(INPUT_FILE_TYPE==FA) 
				Score+=Match_FA_Score;
			else
			{
				float Q_Value=*Qual-QUALITYCONVERSIONFACTOR;//make quality to integer..
				assert(Q_Value>=0);
				float Penalty= -10*log10(Pr(Q_Value));
				Penalty=std::min(QLIMIT_FLOAT,Penalty);
				assert(Penalty<=QLIMIT); 
				Score+= Penalty;
				BQScore+=MATCH_BONUS;
				//QScore+=Penalty;
			}
			SW_Score+=match;
		}
	}
	A.BQScore+=BQScore;
	A.Score+=Score;
	A.QScore+=QScore;
	A.Mismatch+=Mis;
	A.SW_Score+=SW_Score;
}

int Mismatch_Scan(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Ref,char* Read,int Read_Len,int Max_Mis,int Current_Mis)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	int Mis=0;
	int Score=0;
	for(int j=0;j<Read_Len && Mis<=Max_Mis;j++)
	{
		if(*(Ref+j)!=*(Read+j) )//|| (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='T' && (*(rawRead+j)=='c' || *(rawRead+j)=='C' ) ) ||  (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='A' && (*(rawRead+j)=='g' ||*(rawRead+j)=='G' ) ) ) 
		{
			Mis++;
			Score-=mismatch;
		}
		else
		{
			Score+=match;
		}
	}
	return Score;

}

int MatchX(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Ref,int Ref_Len,char* Read,int Read_Len,int Gap,int Max_Mis,Match_Info Pos[],int Current_Mis)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	
	int Least_Mis=INT_MAX,Pos_Ptr=0;
	for(int i=0;(i<=Gap);i++)
	{
		int Mis=Current_Mis;//0;
		int j=0;
		for(;j<Read_Len && Mis<=Max_Mis;j++)
		{
			if(*(Ref+j)!=*(Read+i+j) )//|| (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='T' && (*(rawRead+i+j)=='c' || *(rawRead+i+j)=='C' ) ) ||  (*(Ref+j)==*(rawRef+i+j) && *(rawRef+j)=='A' && (*(rawRead+i+j)=='g' ||*(rawRead+i+j)=='G' ) ) ) 
			{Mis++;}
		}
		if(j==Read_Len && Mis<=Max_Mis)
		{
		       	if(Least_Mis>Mis) 
				Least_Mis=Mis;
			Pos[Pos_Ptr].Length=Read_Len;
			Pos[Pos_Ptr].Mis=Mis;
			Pos[Pos_Ptr++].Loc=i;
		}
		//Read_Len--;
		if(!Read_Len) break;
		assert(Read_Len);
	}
	Pos[Pos_Ptr].Loc=INT_MAX;
	if(Least_Mis<=Max_Mis) return Pos_Ptr;
	else return 0;
}

int Match(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Ref,int Ref_Len,char* Read,int Read_Len,int Gap,int Max_Mis,Match_Info Pos[],int Current_Mis)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	
	int Least_Mis=INT_MAX,Pos_Ptr=0;
	for(int i=0;(i<=Gap);i++)
	{
		int Mis=Current_Mis;//0;
		int j=0;
		for(;j<Read_Len && Mis<=Max_Mis;j++)
		{
			if(*(Ref+j)!=*(Read+i+j))// || (*(Ref+j)==*(rawRef+j) && *(rawRef+j)=='T' && (*(rawRead+i+j)=='c' || *(rawRead+i+j)=='C' ) ) ||  (*(Ref+j)==*(rawRef+i+j) && *(rawRef+j)=='A' && (*(rawRead+i+j)=='g' ||*(rawRead+i+j)=='G' ) ) ) 
			{Mis++;}
		}
		if(j==Read_Len && Mis<=Max_Mis)
		{
		       	if(Least_Mis>Mis) 
				Least_Mis=Mis;
			Pos[Pos_Ptr].Length=Read_Len;
			Pos[Pos_Ptr].Mis=Mis;
			Pos[Pos_Ptr++].Loc=i;
		}
		Read_Len--;
		if(!Read_Len) break;
		assert(Read_Len);
	}
	Pos[Pos_Ptr].Loc=INT_MAX;
	if(Least_Mis<=Max_Mis) return Pos_Ptr;
	else return 0;
}


Gap_Info Get_Max(int Sspan,int Lspan,char Sign,int ReadSpan,READ & RawR,char* rawRef,char* S_Start,char* S_End,char* L_Start,char* L_End,Gap_Info *G_Array,int Mismatch_Limit)
{
/*	READ RawRs;
	if(Sign=='-')
	{
		int i;int Lens=RawR.Real_Len;
		for (i=0;i<=Lens-1;i++){RawRs.Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Tag_Copy[i]];};
	}
	else
	{
		RawRs=RawR;
	}
*/	char* rawRead =RawR.Tag_Copy;//momo
	rawRead+=ReadSpan;
	
	int Cumulative_Mis_Count[L_End-L_Start+1],Cumulative_Mis[L_End-L_Start+1],Mis=0;
	Gap_Info G;
	int Max_Mis_Init=0;//printf("\nErr %s %s %s %s %s %c\n",L_Start,L_End,rawRef,rawRead,S_Start,Sign);
	for(int i=0;L_Start+i<=L_End;i++)
	{
		if(*(S_Start+i)!=*(L_Start+i))//  || (*(L_Start+i)==*(rawRef+i) && *(rawRef+i)=='T' && (*(rawRead+i)=='c' || *(rawRead+i)=='C' ) ) ||  (*(L_Start+i)==*(rawRef+i) && *(rawRef+i)=='A' && (*(rawRead+i)=='g' ||*(rawRead+i)=='G' ) ) ) 
		{
			Mis-=mismatch;
			Max_Mis_Init++;
		}
		else Mis+=match;
		Cumulative_Mis[i]=Mis;
		Cumulative_Mis_Count[i]=Max_Mis_Init;
	}

	int L_Len=(L_End-L_Start)+1;
	int Max_Score=Cumulative_Mis[L_Len-1];Mis=0;int Split_Loc=L_Len;
	int Tot_Score;
	int Max_Mis=0;

	int G_Ptr=0;
	if(Max_Mis_Init<=Mismatch_Limit)
	{
		G_Array[0].Start=Split_Loc;
		G_Array[0].End=L_Len-Split_Loc;
		G_Array[0].Score=Max_Score;
		G_Array[0].Mis= Max_Mis_Init;
		G_Ptr=1;
	}

	for(int i=L_Len-1;i>=0;i--)
	{
		if(*(L_Start+i)!=*(S_End)) 
		{
			Mis-=mismatch;
			Max_Mis++;
			if(Max_Mis>Mismatch_Limit) break;
		}
		else Mis+=match;
		S_End--;
		if(i)
		{
			Tot_Score=Cumulative_Mis[i-1]+Mis;
		}
		else
		{
			Tot_Score=Mis;
		}
		if(Max_Score<=Tot_Score+SCOREDIFF &&  (Max_Mis+(i? Cumulative_Mis_Count[i-1]:0) <=Mismatch_Limit)) 
		{
			if(Max_Score<Tot_Score && !OUTPUTALL)
			{
				G_Ptr=0;
			}

			Max_Score=Tot_Score;
			Split_Loc=i;

			G_Array[G_Ptr].Start=Split_Loc;
			G_Array[G_Ptr].End=L_Len-Split_Loc;
			G_Array[G_Ptr].Score=Max_Score;
			G_Array[G_Ptr++].Mis= Max_Mis+(i? Cumulative_Mis_Count[i-1]:0);
			assert(G_Ptr<=L_Len+1);
		}
	}
	
	G_Array[G_Ptr].Start=INT_MAX;
	G.Start=Split_Loc;
	G.End=L_Len-Split_Loc;
	G.Score=Max_Score;
	G.Mis= (Split_Loc==L_Len) ? Max_Mis_Init:Max_Mis;
	return G;
}

Gap_Info Maximize_Ins(int Sspan,int Lspan,char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Read,char* Ref,char* Ref_End,char* Read_Init,Match_Info Pos[],int Delete,int Offset,int Mismatch_Limit,int Mis_In_Left,int Tot_Read_Len)
{
	int Pos_Ptr=0;
	Gap_Info Final_Gap_Alignment;
	Final_Gap_Alignment.Start=INT_MAX;
	Final_Gap_Alignment.Score=0;
	int Opt_Score=0;
	while(Maximize_Ins < sizeof(Pos)/sizeof(Match_Info) && Pos[Pos_Ptr].Loc!=INT_MAX)
	{
		int Loc=Pos[Pos_Ptr].Loc,Sag=Pos[Pos_Ptr].Length,Mis=Pos[Pos_Ptr++].Mis;

		Gap_Info G[Ref_End-Ref+2];
		Get_Max(Sspan+Loc-1,Lspan,Sign,ReadSpan,RawR,rawRef,Ref,Ref_End,Read,Read_Init+Loc-1,G,Mismatch_Limit);
		int j=0;
		while(G[j].Start!=INT_MAX)
		{
			int Gap_Size=(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End;
			int Matched=Tot_Read_Len-G[j].Mis-Mis-Mis_In_Left-Gap_Size;
			int Final_Score=match*(Matched)-mismatch*(G[j].Mis+Mis+Mis_In_Left)-(gap_open+gap_extension*(Gap_Size-1));
			if(Opt_Score<Final_Score)
			{
				Opt_Score=Final_Score;
				Final_Gap_Alignment=G[j];
				Final_Gap_Alignment.A=G[j].Start+Offset;
				Final_Gap_Alignment.B=(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End;
				Final_Gap_Alignment.Indel_Flag='I';
				Final_Gap_Alignment.C=G[j].End+Sag;
				Final_Gap_Alignment.Mis=G[j].Mis+Mis;
				Final_Gap_Alignment.Score=Final_Score;
			}
			//printf("%dM%dI%dM-%d,%d\n",G[j].Start+Offset,(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
			j++;
		}

	}
	Final_Gap_Alignment.Score=Opt_Score;
	return Final_Gap_Alignment;
}


Gap_Info Maximize_Del(int Sspan,int Lspan,char Sign,int ReadSpan,READ & RawR,char* rawRef,char* Read,char* Ref,char* Ref_End,char* Read_Init,Match_Info Pos[],int Delete,int Offset,int Mismatch_Limit,int Mis_In_Left,int Tot_Read_Len)
{
	int Pos_Ptr=0;
	Gap_Info Gap_Alignment;
	Gap_Alignment.Start=INT_MAX;Gap_Alignment.Score=0;
	int Opt_Score=0;
	while(Pos[Pos_Ptr].Loc!=INT_MAX)
	{
		int Loc=Pos[Pos_Ptr].Loc,Sag=Pos[Pos_Ptr].Length,Mis=Pos[Pos_Ptr++].Mis;

		Gap_Info G[Read_Init+Loc-Read+2];
		Get_Max(Sspan,Lspan+Loc-1,Sign,ReadSpan,RawR,rawRef,Ref,Ref_End+Loc-1,Read,Read_Init,G,Mismatch_Limit);
		int j=0;
		while(G[j].Start!=INT_MAX)
		{
			int Del_Size=(Ref_End+Loc-Ref)-G[j].Start-G[j].End;
			int Matched=Tot_Read_Len-G[j].Mis-Mis-Mis_In_Left-Del_Size;
			int Final_Score=match*(Matched)-mismatch*(G[j].Mis+Mis+Mis_In_Left)-(gap_open+gap_extension*(Del_Size-1));
			if(Opt_Score<Final_Score)
			{
				Opt_Score=Final_Score;
				Gap_Alignment=G[j];
				Gap_Alignment.A=G[j].Start+Offset;
				Gap_Alignment.B=Del_Size;
				Gap_Alignment.Indel_Flag='D';
				Gap_Alignment.C=G[j].End+Sag;
				Gap_Alignment.Mis=G[j].Mis+Mis;
				Gap_Alignment.Score=Final_Score;
			}
			//printf("%dM%dD%dM-%d,%d\n",G[j].Start+Offset,Del_Size,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
			j++;
		}
	}
	return Gap_Alignment;
}

int Insert_In_First_Half(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* D,char* S,int Delete,int Insert,int Max_Mis,int Offset,int Mis_In_Left,int Tot_Read_Len,Gap_Info & Gap_Alignment)
{
	int LS=strlen(S);
	int Backtrack=LS/2;
	Gap_Info G;

	int Gap=Insert;
	if(Gap >LS/2) 
	{
		Gap=LS/2;//Insert=Gap-Delete;
	}
	Match_Info Pos[Gap+1],Pos2[Gap+1],Pos1[Gap+1];

	int S1=Match(Sign,LS/2+ReadSpan,RawR,rawRef+LS/2,D+LS/2,strlen(D),S+Backtrack,(LS-Backtrack)/*LS/2+Delete)-Length*/,Gap,Max_Mis,Pos,Mis_In_Left);//Indel in first third
		
	assert(Max_Mis-Mis_In_Left>=0);
	if(S1)
	{
		//return Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,LS/2,Delete,Offset);
		G=Maximize_Ins(Backtrack,LS/2-1,Sign,ReadSpan,RawR,rawRef,S,D,D+LS/2-1,S+Backtrack,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
		if(G.Score>Gap_Alignment.Score)
		{
			Gap_Alignment=G;
		}
	}
	if(Insert>=LS/2)
	{
		Pos[0].Length=0;
		Pos[0].Mis=0;
		Pos[0].Loc=0;
		Pos[1].Loc=INT_MAX;
		for(int i=0;i<LS/2;i++)
		{
			G=Maximize_Ins(LS,i+1,Sign,ReadSpan,RawR,rawRef,S,D,D+i+1,S+LS,Pos,0,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
			if(G.Score>Gap_Alignment.Score)
			{
				Gap_Alignment=G;
			}
		}
	}
	return 1;
};

int Delete_In_First_Half(char Sign,int ReadSpan,READ & RawR,char* rawRef,char* D,char* S,int Delete,int Insert,int Max_Mis,int Offset,int Mis_In_Left,int Tot_Read_Len,Gap_Info & Gap_Alignment)
{
	int LS=strlen(S);
	int Backtrack=LS/2;

	int Gap=Delete;
	Match_Info Pos[Gap+1],Pos2[Gap+1],Pos1[Gap+1];

	int S1=MatchX(Sign,LS/2+ReadSpan,RawR,rawRef+LS/2,D+LS/2,(LS-LS/2)/*LS/2+Delete)-Length*/,S+LS/2,strlen(S),Gap,Max_Mis,Pos,Mis_In_Left);//Indel in first third
	assert(Max_Mis-Mis_In_Left>=0);		
	if(S1)
	{
		//return Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,LS/2,Delete,Offset);
		Gap_Info G=Maximize_Del(LS/2-1,LS/2,Sign,ReadSpan,RawR,rawRef,S,D,D+LS/2,S+LS/2-1,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);

		if(Gap_Alignment.Score<G.Score)
		{
			Gap_Alignment=G;
		}
	}
	/*if(Delete>=LS/2)
	{
		Pos[0].Length=0;
		Pos[0].Mis=0;
		Pos[0].Loc=0;
		Pos[1].Loc=INT_MAX;
		for(int i=0;i<LS/2;i++)
		{
			Maximize_Del(S,D,D+i+1,S+LS,Pos,0,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
		}
	}*/
	return 1;
};

int Ham(char* D,char* S,int Len)
{
	int Mis=0;
	for(int i=0;i<Len;i++)
	{
		if(*(S+i)!=*(D+i)) Mis++;
	}
	return Mis;
};

extern int INSERT;
extern int DELETE;
int Max_Mis=2;
int Max_Score;
//const int OFF=75;

Alignment Fast_SW(char *Ref,char *Read,char *Quality,int OFF,char Sign,char* rawRef,READ & RawR)
{
	if(OFF==INT_MAX) OFF=SEEDSIZE;
	int Tot_Read_Len_Org=strlen(Read);
	int Tot_Read_Len=strlen(Read);
	//char *D=Ref+Tot_Read_Len/2,*S=Read+Tot_Read_Len/2,*Q=Quality+Tot_Read_Len/2;
	char *D=Ref+OFF,*S=Read+OFF,*Q=Quality+OFF;
	int Offset=0;
	int Read_Len;

	int Mis_In_LH=0;
	int Score_In_LH=0;
	//int Half_Length=Tot_Read_Len-Tot_Read_Len/2;
	int Half_Length=Tot_Read_Len-OFF;
	Gap_Info Gap_Alignment;Gap_Alignment.Score= -1;
	int TScore=Mismatch_Scan(Sign,OFF,RawR,rawRef+OFF,D,S,Half_Length,100,0);
	if(TScore >0)
	{
		Gap_Alignment.Score=TScore;
		Gap_Alignment.A=Tot_Read_Len_Org;Gap_Alignment.B=Gap_Alignment.C=0;
	}
	do
	{
		//Max_Score=Insert_In_First_Half(D,S,Delete,Insert,Max_Mis,Offset,Mis_In_LH,Tot_Read_Len);
		Max_Score=Delete_In_First_Half(Sign,OFF,RawR,rawRef+OFF,D,S,DELETE,INSERT,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
		Read_Len=strlen(S);
		Mis_In_LH+=Ham(D,S,Read_Len/2);
		if(Mis_In_LH>Max_Mis) break;
		//Score_In_LH=Mis_In_LH*MISMATCH+(Tot_Read_Len_LH-Max_Mis)*MATCH;
		D=D+Read_Len/2;S=S+Read_Len/2;
		Offset+=Read_Len/2;
	}
	while ((Read_Len-Read_Len/2)>1);

	//D=Ref+Tot_Read_Len/2;S=Read+Tot_Read_Len/2;Q=Quality+Tot_Read_Len/2;
	D=Ref+OFF,S=Read+OFF,Q=Quality+OFF;
	Offset=0;
	Mis_In_LH=0;
	Score_In_LH=0;
	do
	{
		Max_Score=Insert_In_First_Half(Sign,OFF,RawR,rawRef+OFF,D,S,DELETE,INSERT,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
		Read_Len=strlen(S);
		Mis_In_LH+=Ham(D,S,Read_Len/2);
		if(Mis_In_LH>Max_Mis) break;
		//Score_In_LH=Mis_In_LH*MISMATCH+(Tot_Read_Len_LH-Max_Mis)*MATCH;
		D=D+Read_Len/2;S=S+Read_Len/2;
		Offset+=Read_Len/2;
	}
	while ((Read_Len-Read_Len/2)>1);
	Alignment A;
	A.Loc=0;
	if(Gap_Alignment.Score>0 && Gap_Alignment.Score!=INT_MAX)
	{
		A.SW_Score=0;A.Mismatch=0;
		A.Indel=Gap_Alignment.B;A.BQScore=A.Score=A.QScore=A.Clip_H=A.Clip_T=0;
		if(Gap_Alignment.B)//Indel present
		{
			int Skip=Tot_Read_Len-Half_Length+Gap_Alignment.A;//skip from begin..
			Mismatch_Scan_With_Score(Sign,0,RawR,Ref,Read,Quality,Skip,100,0,A,rawRef);//left half..
			if(Gap_Alignment.Indel_Flag=='I')
			{
				Mismatch_Scan_With_Score(Sign,Skip+Gap_Alignment.B,RawR,Ref+Skip,Read+Skip+Gap_Alignment.B,Quality+Skip+Gap_Alignment.B,Tot_Read_Len-Skip-Gap_Alignment.B,100,0,A,rawRef);
				A.Loc+=Gap_Alignment.B;
			}
			else
			{
				assert(Gap_Alignment.Indel_Flag=='D');
				Mismatch_Scan_With_Score(Sign,Skip,RawR,Ref+Skip+Gap_Alignment.B,Read+Skip,Quality+Skip,Tot_Read_Len-Skip,100,0,A,rawRef);
				A.Loc-=Gap_Alignment.B;
			}
			A.BQScore-=(BOPEN+BEXT*(Gap_Alignment.B));
			A.SW_Score-=(gap_open+Gap_Alignment.B*gap_extension);
			A.Score+=(gap_openP+Gap_Alignment.B*gap_extensionP);
			if(Sign=='+')
			{
				sprintf(A.Cigar,"%dM%d%c%dM",Skip,Gap_Alignment.B,Gap_Alignment.Indel_Flag,Gap_Alignment.C);
			}
			else
			{
				sprintf(A.Cigar,"%dM%d%c%dM",Gap_Alignment.C,Gap_Alignment.B,Gap_Alignment.Indel_Flag,Skip);
			}
		}
		else//Mismatch hit..
		{
			Mismatch_Scan_With_Score(Sign,0,RawR,Ref,Read,Quality,Tot_Read_Len,100,0,A,rawRef);
			sprintf(A.Cigar,"%dM",Gap_Alignment.A);
		}
		A.Score= -A.Score;
	}
	else
	{
		assert(Gap_Alignment.Score<=0);
		A.Score=INT_MAX;
		A.SW_Score=0;
	}
	return A;
}

Alignment Fast_SWX(char Sign,READ & RawR,char* rawRef,char *Ref,char *Read,int OFF)
{
	if(OFF==INT_MAX) OFF=SEEDSIZE;
	int Tot_Read_Len_Org=strlen(Read);
	int Tot_Read_Len=strlen(Read);
	//char *D=Ref+Tot_Read_Len/2,*S=Read+Tot_Read_Len/2,*Q=Quality+Tot_Read_Len/2;
	char *D=Ref+OFF,*S=Read+OFF,*rawRefD=rawRef+OFF;
	int Offset=0;
	int Read_Len;

	int Mis_In_LH=0;
	int Score_In_LH=0;
	//int Half_Length=Tot_Read_Len-Tot_Read_Len/2;
	int Half_Length=Tot_Read_Len-OFF;
	Gap_Info Gap_Alignment;Gap_Alignment.Score= -1;
	int TScore=Mismatch_Scan(Sign,OFF,RawR,rawRefD,D,S,Half_Length,100,0);
	if(TScore >0)
	{
		Gap_Alignment.Score=TScore;
		Gap_Alignment.A=Tot_Read_Len_Org;Gap_Alignment.B=Gap_Alignment.C=0;
	}
	do
	{
		//Max_Score=Insert_In_First_Half(D,S,Delete,Insert,Max_Mis,Offset,Mis_In_LH,Tot_Read_Len);
		Max_Score=Delete_In_First_Half(Sign,OFF,RawR,rawRefD,D,S,DELETE,INSERT,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
		Read_Len=strlen(S);
		Mis_In_LH+=Ham(D,S,Read_Len/2);
		if(Mis_In_LH>Max_Mis) break;
		//Score_In_LH=Mis_In_LH*MISMATCH+(Tot_Read_Len_LH-Max_Mis)*MATCH;
		D=D+Read_Len/2;S=S+Read_Len/2;
		rawRefD=rawRefD+Read_Len/2;OFF+=Read_Len/2;//moxian
		Offset+=Read_Len/2;
	}
	while ((Read_Len-Read_Len/2)>1);

	D=Ref+OFF,S=Read+OFF;
	Offset=0;
	Mis_In_LH=0;
	Score_In_LH=0;
	do
	{
		Max_Score=Insert_In_First_Half(Sign,OFF,RawR,rawRefD,D,S,DELETE,INSERT,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
		Read_Len=strlen(S);
		Mis_In_LH+=Ham(D,S,Read_Len/2);
		if(Mis_In_LH>Max_Mis) break;
		//Score_In_LH=Mis_In_LH*MISMATCH+(Tot_Read_Len_LH-Max_Mis)*MATCH;
		D=D+Read_Len/2;S=S+Read_Len/2;
		rawRefD=rawRefD+Read_Len/2;OFF+=Read_Len/2;//moxian
		Offset+=Read_Len/2;
	}
	while ((Read_Len-Read_Len/2)>1);
	Alignment A;
	if(Gap_Alignment.Score>0 && Gap_Alignment.Score!=INT_MAX)
	{
		A.SW_Score=0;A.Mismatch=0;
		A.Indel=Gap_Alignment.B;A.BQScore=A.Score=A.QScore=A.Clip_H=A.Clip_T=0;
		if(Gap_Alignment.B)//Indel present
		{
			int Skip=Tot_Read_Len-Half_Length+Gap_Alignment.A;//skip from begin..
			Mismatch_Scan_With_ScoreX(Sign,0,RawR,rawRef,Ref,Read,Skip,100,0,A);//left half..
			if(Gap_Alignment.Indel_Flag=='I')
			{
				Mismatch_Scan_With_ScoreX(Sign,Skip+Gap_Alignment.B,RawR,rawRef+Skip,Ref+Skip,Read+Skip+Gap_Alignment.B,Tot_Read_Len-Skip-Gap_Alignment.B,100,0,A);
			}
			else
			{
				assert(Gap_Alignment.Indel_Flag=='D');
				//Mismatch_Scan_With_Score(Ref+Skip+Gap_Alignment.B,Read+Skip,Quality,Tot_Read_Len-Skip+Gap_Alignment.B,100,0,A);
				Mismatch_Scan_With_ScoreX(Sign,Skip,RawR,rawRef+Skip+Gap_Alignment.B,Ref+Skip+Gap_Alignment.B,Read+Skip,Tot_Read_Len-Skip,100,0,A);
			}
			A.BQScore-=(BOPEN+BEXT*(Gap_Alignment.B));
			A.SW_Score-=(gap_open+Gap_Alignment.B*gap_extension);
			A.Score+=(gap_openP+Gap_Alignment.B*gap_extensionP);
			sprintf(A.Cigar,"%dM%d%c%dM",Skip,Gap_Alignment.B,Gap_Alignment.Indel_Flag,Gap_Alignment.C);
		}
		else//Mismatch hit..
		{
			Mismatch_Scan_With_ScoreX(Sign,0,RawR,rawRef,Ref,Read,Tot_Read_Len,100,0,A);
			sprintf(A.Cigar,"%dM",Gap_Alignment.A);
		}
		A.Score= -A.Score;
	}
	else
	{
		assert(Gap_Alignment.Score<=0);
		A.Score=INT_MAX;
	}
	return A;
}
