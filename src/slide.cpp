#include <cstdio>
#include <string>
#include "string.h"
#include "limits.h"
#include "assert.h"
#include "fstream"
#include "iostream"
#include "math.h"
#include <cstdlib>
#include <map>
using namespace std;
const bool OUTPUTALL=true;
const int MATCH=2;
const int MISMATCH=-2;
const int OPEN=-4;
const int EXTEND=-1;
const int SCOREDIFF=8*MATCH;
const int MAX_SIGLEN=10;
const int gap_open = 6; //原来是3改为6
const int gap_extension = 1;
int MX=6,MN=2,BOPEN=6,BEXT=1,MATCH_BONUS=0;
int gap_openP=40,gap_extensionP=6; 
float Match_FA_Score;
float Mis_FA_Score; 
struct Alignment
{
	unsigned Loc;
	char Sign;
	//unsigned char Bad_Loc[SW_MAX_MIS_TO_STORE];
	//char Bad_Char[10];
	//char Edit_Distance;
	int Indel;
	int Mismatch;
	int Score;
	int QScore;
	int SW_Score;
	int BQScore;
	int Sub_Opt_Score;
	char Realigned;
	char Cigar[MAX_SIGLEN+1];
	int Clip_H,Clip_T;
};


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
const int FA=1;
const float QLIMIT_FLOAT=30;
const float QLIMIT=QLIMIT_FLOAT;
const int QUALITYCONVERSIONFACTOR=30;
const int INPUT_FILE_TYPE=0; 
const int POWLIMIT=300;
float POW10[POWLIMIT];

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

void Mismatch_Scan_With_Score(char* Ref,char* Read,char* Qual,int Read_Len,int Max_Mis,int Current_Mis,Alignment & A)
{
	int Mis=0;
	int SW_Score=0;
	float Score=0,QScore=0,BQScore=0;
	float TScore=0,TQScore=0,TBQScore=0;
	for(int j=0;j<Read_Len && Mis<=Max_Mis;j++,Qual++)
	{
		if(*(Ref+j)!=*(Read+j)) 
		{
			Mis++;
			SW_Score+=MISMATCH;
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
			SW_Score+=MATCH;
		}
	}
	A.BQScore+=BQScore;
	A.Score+=Score;
	A.Mismatch+=Mis;
	A.SW_Score+=SW_Score;
}

int Mismatch_Scan(char* Ref,char* Read,int Read_Len,int Max_Mis,int Current_Mis)
{
	int Mis=0;
	int Score=0;
	for(int j=0;j<Read_Len && Mis<=Max_Mis;j++)
	{
		if(*(Ref+j)!=*(Read+j)) 
		{
			Mis++;
			Score+=MISMATCH;
		}
		else
		{
			Score+=MATCH;
		}
	}
	return Score;

}

int MatchX(char* Ref,int Ref_Len,char* Read,int Read_Len,int Gap,int Max_Mis,Match_Info Pos[],int Current_Mis)
{
	int Least_Mis=INT_MAX,Pos_Ptr=0;
	for(int i=0;(i<=Gap);i++)
	{
		int Mis=Current_Mis;//0;
		int j=0;
		for(;j<Read_Len && Mis<=Max_Mis;j++)
		{
			if(*(Ref+j)!=*(Read+i+j)) Mis++;
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

int Match(char* Ref,int Ref_Len,char* Read,int Read_Len,int Gap,int Max_Mis,Match_Info Pos[],int Current_Mis)
{
	int Least_Mis=INT_MAX,Pos_Ptr=0;
	for(int i=0;(i<=Gap);i++)
	{
		int Mis=Current_Mis;//0;
		int j=0;
		for(;j<Read_Len && Mis<=Max_Mis;j++)
		{
			if(*(Ref+j)!=*(Read+i+j)) Mis++;
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


Gap_Info Get_Max(char* S_Start,char* S_End,char* L_Start,char* L_End,Gap_Info *G_Array,int Mismatch_Limit)
{

	int Cumulative_Mis_Count[L_End-L_Start+1],Cumulative_Mis[L_End-L_Start+1],Mis=0;
	Gap_Info G;
	int Max_Mis_Init=0;
	for(int i=0;L_Start+i<=L_End;i++)
	{
		if(*(S_Start+i)!=*(L_Start+i)) 
		{
			Mis+=MISMATCH;
			Max_Mis_Init++;
		}
		else Mis+=MATCH;
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
			Mis+=MISMATCH;
			Max_Mis++;
			if(Max_Mis>Mismatch_Limit) break;
		}
		else Mis+=MATCH;
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

Gap_Info Maximize_Ins(char* Read,char* Ref,char* Ref_End,char* Read_Init,Match_Info Pos[],int Delete,int Offset,int Mismatch_Limit,int Mis_In_Left,int Tot_Read_Len)
{
	int Score,Init_Score=0;
	int Pos_Ptr=0;
	Gap_Info Final_Gap_Alignment;
	Final_Gap_Alignment.Start=INT_MAX;
	Final_Gap_Alignment.Score=0;
	int Opt_Score=0;
	while(Pos[Pos_Ptr].Loc!=INT_MAX)
	{
		int Loc=Pos[Pos_Ptr].Loc,Sag=Pos[Pos_Ptr].Length,Mis=Pos[Pos_Ptr++].Mis;
		Score=Init_Score+MISMATCH*Mis;

		Gap_Info G[Ref_End-Ref+2];
		Get_Max(Read,Read_Init+Loc-1,Ref,Ref_End,G,Mismatch_Limit);
		int j=0;
		while(G[j].Start!=INT_MAX)
		{
			int Gap_Size=(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End;
			int Matched=Tot_Read_Len-G[j].Mis-Mis-Gap_Size;
			int Final_Score=MATCH*(Matched)+MISMATCH*(G[j].Mis+Mis)+(OPEN+EXTEND*(Gap_Size-1));
			if(Opt_Score<Final_Score)
			{
				Opt_Score=Final_Score;
				Final_Gap_Alignment=G[j];
				Final_Gap_Alignment.A=G[j].Start+Offset;
				Final_Gap_Alignment.B=(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End;
				Final_Gap_Alignment.Indel_Flag='I';
				Final_Gap_Alignment.C=G[j].End+Sag;
				Final_Gap_Alignment.Mis=G[j].Mis+Mis;
			}
			//printf("%dM%dI%dM-%d,%d\n",G[j].Start+Offset,(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
			j++;
		}

	}
	Final_Gap_Alignment.Score=Opt_Score;
	return Final_Gap_Alignment;
}


Gap_Info Maximize_Del(char* Read,char* Ref,char* Ref_End,char* Read_Init,Match_Info Pos[],int Delete,int Offset,int Mismatch_Limit,int Mis_In_Left,int Tot_Read_Len)
{
	int Score,Init_Score=0;
	int Pos_Ptr=0;
	Gap_Info Gap_Alignment;
	Gap_Alignment.Start=INT_MAX;Gap_Alignment.Score=0;
	int Opt_Score=0;
	while(Pos[Pos_Ptr].Loc!=INT_MAX)
	{
		int Loc=Pos[Pos_Ptr].Loc,Sag=Pos[Pos_Ptr].Length,Mis=Pos[Pos_Ptr++].Mis;
		Score=Init_Score+MISMATCH*Mis;

		Gap_Info G[Read_Init+Loc-Read+2];
		Get_Max(Ref,Ref_End+Loc-1,Read,Read_Init,G,Mismatch_Limit);
		int j=0;
		while(G[j].Start!=INT_MAX)
		{
			int Del_Size=(Ref_End+Loc-Ref)-G[j].Start-G[j].End;
			int Matched=Tot_Read_Len-G[j].Mis-Mis-Del_Size;
			int Final_Score=MATCH*(Matched)+MISMATCH*(G[j].Mis+Mis)+(OPEN+EXTEND*(Del_Size-1));
			if(Opt_Score<Final_Score)
			{
				Opt_Score=Final_Score;
				Gap_Alignment=G[j];
				Gap_Alignment.A=G[j].Start+Offset;
				Gap_Alignment.B=Del_Size;
				Gap_Alignment.Indel_Flag='D';
				Gap_Alignment.C=G[j].End+Sag;
				Gap_Alignment.Mis=G[j].Mis+Mis;
			}
			//printf("%dM%dD%dM-%d,%d\n",G[j].Start+Offset,Del_Size,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
			j++;
		}
	}
	return Gap_Alignment;
}

/*int Maximize(char* Read,char* Ref,char* Ref_End,char* Read_Init,Match_Info Pos[],int Delete,int Offset,int Mismatch_Limit,int Mis_In_Left,int Tot_Read_Len)
{
	int Score,Init_Score=0;
	int Pos_Ptr=0;
	while(Pos[Pos_Ptr].Loc!=INT_MAX)
	{
		int Loc=Pos[Pos_Ptr].Loc,Sag=Pos[Pos_Ptr].Length,Mis=Pos[Pos_Ptr++].Mis;
		Score=Init_Score+MISMATCH*Mis;
		if(Loc<Delete)
		{
			Gap_Info G[Read_Init+Loc-Read];
			Get_Max(Ref,Ref_End,Read,Read_Init+Loc-1,G,Mismatch_Limit);
			int j=0;
			while(G[j].Start!=INT_MAX)
			{
				int Final_Score=MATCH*(Tot_Read_Len-G[j].Mis-Mis)+MISMATCH*(G[j].Mis+Mis)+(OPEN+EXTEND*(Loc)+G[j].Score);
				//int Final_Score=Score+(OPEN+EXTEND*(Loc)+G[j].Score);
				printf("%dM%dD%dM-%d,%d\n",G[j].Start+Offset,(Ref_End-Ref+1)-G[j].Start-G[j].End,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
				j++;
			}
		}
		else
		{
			Gap_Info G[Ref_End-Ref+2];
			Get_Max(Read,Read_Init+Loc-1,Ref,Ref_End,G,Mismatch_Limit);
			int j=0;
			while(G[j].Start!=INT_MAX)
			{
				int Gap_Size=(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End;
				int Matched=Tot_Read_Len-G[j].Mis-Mis-Gap_Size;
				int Final_Score=MATCH*(Matched)+MISMATCH*(G[j].Mis+Mis)+(OPEN+EXTEND*(Gap_Size-1));
				printf("%dM%dI%dM-%d,%d\n",G[j].Start+Offset,(Read_Init+Loc-1-Read)+1-G[j].Start-G[j].End,G[j].End+Sag,G[j].Mis+Mis,Final_Score);
				j++;
			}
		}

	}
}*/

int Insert_In_First_Half(char* D,char* S,int Delete,int Insert,int Max_Mis,int Offset,int Mis_In_Left,int Tot_Read_Len,Gap_Info & Gap_Alignment)
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

	int S1=Match(D+LS/2,strlen(D),S+Backtrack,(LS-Backtrack)/*LS/2+Delete)-Length*/,Gap,Max_Mis,Pos,Mis_In_Left);//Indel in first third
	assert(Max_Mis-Mis_In_Left>=0);
	if(S1)
	{
		//return Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,LS/2,Delete,Offset);
		G=Maximize_Ins(S,D,D+LS/2-1,S+Backtrack,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
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
			G=Maximize_Ins(S,D,D+i+1,S+LS,Pos,0,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
			if(G.Score>Gap_Alignment.Score)
			{
				Gap_Alignment=G;
			}
		}
	}
	return 1;
};

int Delete_In_First_Half(char* D,char* S,int Delete,int Insert,int Max_Mis,int Offset,int Mis_In_Left,int Tot_Read_Len,Gap_Info & Gap_Alignment)
{
	int LS=strlen(S);
	int Backtrack=LS/2;

	int Gap=Delete;
	Match_Info Pos[Gap+1],Pos2[Gap+1],Pos1[Gap+1];

	int S1=MatchX(S+LS/2,strlen(S),D+LS/2,(LS-LS/2)/*LS/2+Delete)-Length*/,Gap,Max_Mis,Pos,Mis_In_Left);//Indel in first third
	assert(Max_Mis-Mis_In_Left>=0);
	if(S1)
	{
		//return Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,LS/2,Delete,Offset);
		Gap_Info G=Maximize_Del(S,D,D+LS/2,S+LS/2-1,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
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
/*int Indel_In_First_Half(char* D,char* S,int Delete,int Insert,int Max_Mis,int Offset,int Mis_In_Left,int Tot_Read_Len)
{
	int LS=strlen(S);
	int Backtrack;
	if(Delete>LS/2)
	{
		Delete=LS/2;
	}
	if(Delete<=LS/2)
		Backtrack=LS/2-Delete;

	int Gap=Insert+Delete;
	if(Gap >LS/2) 
	{
		Gap=LS;//Insert=Gap-Delete;
	}
	Match_Info Pos[Gap+1],Pos2[Gap+1],Pos1[Gap+1];

	int S1=Match(D+LS/2,strlen(D),S+Backtrack,(LS-Backtrack),Gap,Max_Mis,Pos,Mis_In_Left);//Indel in first third
	if(Delete==LS/2)
	{
		Pos[S1].Length=0;
		Pos[S1].Mis=Mis_In_Left;
		Pos[S1++].Loc=Gap;
		Pos[S1].Loc=INT_MAX;
	}
	assert(Max_Mis-Mis_In_Left>=0);
	if(S1)
	{
		//return Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,LS/2,Delete,Offset);
		Maximize(S,D,D+LS/2-1,S+Backtrack,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
	}
	if(Insert>=LS/2)
	{
		Pos[0].Length=0;
		Pos[0].Mis=0;
		Pos[0].Loc=0;
		Pos[1].Loc=INT_MAX;
		for(int i=0;i<LS/2;i++)
		{
			Maximize(S,D,D+i+1,S+LS,Pos,0,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
		}
	}
	if(Delete>=LS/2)
	{
		Pos[0].Length=0;
		Pos[0].Mis=0;
		Pos[0].Loc=0;
		Pos[1].Loc=INT_MAX;
		for(int i=0;i<LS;i++)
		{
			Maximize(S,D,D+i+LS+1,S+LS,Pos,Delete,Offset,Max_Mis-Mis_In_Left,Mis_In_Left,Tot_Read_Len);
		}
	}
	return 1;
};*/

int Ham(char* D,char* S,int Len)
{
	int Mis=0;
	for(int i=0;i<Len;i++)
	{
		if(*(S+i)!=*(D+i)) Mis++;
	}
	return Mis;
};

int Insert=5;
int Delete=5;
int Max_Mis=0;
int Max_Score;

Alignment Fast_SW(char *Ref,char *Read,char *Quality);

int main(int argc,char* argv[])
{

	char Des[100];//     "ACGTTTACGTTTACGTTTACGTACGTACGTTTACGTTTAAAAAAAAAAAA";
	char Source[100];
	char Quality[100];

	int LS;//strlen(S);


	if(argc <3)
	{
		cout <<argv[0]<<" Genome.fa Reads.fa insert_size delete_size Mis_allowed\n";
		exit(100);
	}
	if(argc >3)
	{
		Insert=atoi(argv[3]);
		Delete=atoi(argv[4]);
		Max_Mis=atoi(argv[5]);
	}

	ifstream I,GEN;
	GEN.open(argv[1]);
	I.open(argv[2]);

	int Read_No=0;
	string line;
	getline (GEN,line);
	getline (GEN,line);
	strcpy(Des,line.c_str());

	Build_Pow10();
	int DEFAULTFASTAQUAL=35;//40;
	Match_FA_Score= -10*log10(Pr(DEFAULTFASTAQUAL));
	Mis_FA_Score  = -10*log10((1-Pr(DEFAULTFASTAQUAL))/3);

	while ( I.good() )
	{
		getline (I,line);
		cout<< line<<endl;
		getline (I,line);
		strcpy(Source,line.c_str());
		getline (I,line);
		strcpy(Quality,line.c_str());
		Alignment A=Fast_SW(Des,Source,Quality);
		Read_No++;
		printf("%s-%d-%d-%d\n",A.Cigar,A.SW_Score,A.Mismatch,A.Score);

	}
}

const int OFF=10;
Alignment Fast_SW(char *Ref,char *Read,char *Quality)
{
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
	int TScore=Mismatch_Scan(D,S,Half_Length,100,0);
	if(TScore >0)
	{
		Gap_Alignment.Score=TScore;
		Gap_Alignment.A=Tot_Read_Len_Org;Gap_Alignment.B=Gap_Alignment.C=0;
	}
	do
	{
		//Max_Score=Insert_In_First_Half(D,S,Delete,Insert,Max_Mis,Offset,Mis_In_LH,Tot_Read_Len);
		Max_Score=Delete_In_First_Half(D,S,Delete,Insert,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
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
		Max_Score=Insert_In_First_Half(D,S,Delete,Insert,Max_Mis,Offset,Mis_In_LH,Half_Length,Gap_Alignment);
		Read_Len=strlen(S);
		Mis_In_LH+=Ham(D,S,Read_Len/2);
		if(Mis_In_LH>Max_Mis) break;
		//Score_In_LH=Mis_In_LH*MISMATCH+(Tot_Read_Len_LH-Max_Mis)*MATCH;
		D=D+Read_Len/2;S=S+Read_Len/2;
		Offset+=Read_Len/2;
	}
	while ((Read_Len-Read_Len/2)>1);
	assert(Gap_Alignment.Score!= -1);
	Alignment A;
	if(Gap_Alignment.Score!=INT_MAX)
	{
		A.SW_Score=0;A.Mismatch=0;
		A.Indel=Gap_Alignment.B;A.BQScore=A.Score=A.Clip_H=A.Clip_T=0;
		if(Gap_Alignment.B)//Indel present
		{
			int Skip=Tot_Read_Len-Half_Length+Gap_Alignment.A;//skip from begin..
			Mismatch_Scan_With_Score(Ref,Read,Quality,Skip,100,0,A);//left half..
			if(Gap_Alignment.Indel_Flag=='I')
			{
				Mismatch_Scan_With_Score(Ref+Skip,Read+Skip+Gap_Alignment.B,Quality,Tot_Read_Len-Skip-Gap_Alignment.B,100,0,A);
			}
			else
			{
				assert(Gap_Alignment.Indel_Flag=='D');
				Mismatch_Scan_With_Score(Ref+Skip+Gap_Alignment.B,Read+Skip,Quality,Tot_Read_Len-Skip+Gap_Alignment.B,100,0,A);
			}
			A.BQScore+=(OPEN+EXTEND*(Gap_Alignment.B));
			A.SW_Score-=(gap_open+Gap_Alignment.B*gap_extension);
			A.Score+=(gap_openP+Gap_Alignment.B*gap_extensionP);
			sprintf(A.Cigar,"%dM%d%c%dM",Gap_Alignment.A+Tot_Read_Len/2,Gap_Alignment.B,Gap_Alignment.Indel_Flag,Gap_Alignment.C);
		}
		else//Mismatch hit..
		{
			Mismatch_Scan_With_Score(Ref,Read,Quality,Tot_Read_Len,100,0,A);
			sprintf(A.Cigar,"%dM",Gap_Alignment.A);
		}
	}
	return A;
}
