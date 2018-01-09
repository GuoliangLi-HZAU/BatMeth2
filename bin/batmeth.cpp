#include "config.h"
#include <string>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <signal.h>
//#include <dvec.h>
#include "zlib.h"
#include <queue>
#include <ctype.h>
#include "batlib.h"
#include <pthread.h>
#include "map.h"
#include "commandline.h"
#include <string>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <omp.h>
//using namespace std;
using namespace std;
const int FILEBUFSIZE=10000000;//Space for batman to write records.. 600000
const int MAX_REC_SIZE = 2000;//maximum size a record of batman will occupy..
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
int MAXHITS = 195;
int maxhits=150;
//int muchHits=8000;
int allHits=0;
int LOOKUPSIZE;
int MAX_MISMATCHES=2;
bool ORIENTATION=false; //scan all 4 possible orientation of a read by default (Cokus protocol)
int THREAD = 4;
int TRIM_LENGTH = 0;
unsigned MAX_READS=0;
char Char_To_C[255];
std::map <unsigned, Ann_Info> Annotations;
std::map <string,int> String_Hash;
BWT *fwfmiCT,*revfmiCT,*fwfmiGA,*revfmiGA;
MEMLOOK MLookCT,MLookGA;
INFILE F;
FILE* TEMP_FQ;
inline int getMismatch(char* s2);
inline int getNtab(char* s2);
inline void Print_List(bool & Multipe,int & minTrueMis,bool & Header_Printed,int cntHit,int firstmismatch,int lastMismatch,std::string* tmpList,char* Header,FILE *Output_File);
void *Do_Meth( void *T);
bool Get_Records(char* File_Buffer, std::string*  tmpList,int & firstmismatch,int & lastMismatch,int & cntHit,char Count);
void Build_Names(const char* Genome_Name,FMFILES & FM);//build FM index file names...
void Parse_Command_line(int argc, char* argv[],char* & GENOME,char* & CT_GENOME,char* & GA_GENOME,char* & INPUTFILE,char* & OUTPUTFILE,int & MAXHITS, int & MAXMISMATCH,unsigned & MAX_READS);
inline void ReplaceCtoT(char* Read);
inline void ReplaceGtoA(char* Read);
//unsigned Location_ArrayCT[80],Location_ArrayGA[80];
//Offset_Record Genome_OffsetsCT[80],Genome_OffsetsGA[80];
Offset_Record Genome_Offsets[80];
int Genome_CountCT;
int Genome_CountGA;
//LEN L;
char* OUTPUTFILE;
int Genome_CountX=0;
void SolexaMapping(int argc, char* argv[]);
void Launch_Threads(int NTHREAD, void* (*Map_t)(void*),Thread_Arg T);
Gene_Hash Genome_List[80];
//bool Multipe=false;
int minMis=200;

int main(int argc, char* argv[])
{	
	//printf("\nBatMeth: Mapper for Bisulfite-seq Reads v1.04b\n");
	//printf("batmeth -g GENOME -i METH -n <number of mismatches> -o <Output_File> -p <#thread>\n\n");

	
	SolexaMapping(argc, argv);

}

void SolexaMapping(int argc, char* argv[]) 
{
	char *CT_GENOME,*GA_GENOME,*GENOME;
	char* INPUTFILE;
	Parse_Command_line(argc, argv,GENOME,CT_GENOME,GA_GENOME,INPUTFILE,OUTPUTFILE,MAXHITS, MAX_MISMATCHES,MAX_READS);
	Char_To_C['A']='T';Char_To_C['a']='t';
	Char_To_C['C']='G';Char_To_C['c']='g';
	Char_To_C['G']='C';Char_To_C['g']='c';
	Char_To_C['T']='A';Char_To_C['t']='a';
	Char_To_C['N']='N';Char_To_C['n']='n';

	char Header[MAX_REC_SIZE];
//{-------------------------------- FM File Stuff ---------------------------------------------------
	FMFILES CT,GA;
	MMPool *mmPool;
//	L.IGNOREHEAD=0;

	std::string CTS=GENOME,GAS=GENOME;
	if (GENOME)
	{
		CTS+="-CtoT";GAS+="-GtoA";
		CT_GENOME=(char*)CTS.c_str();GA_GENOME=(char*)GAS.c_str();
	}

	Build_Names(CT_GENOME,CT);Build_Names(GA_GENOME,GA);
	//batmeth2 
	string G=GENOME;G+=".bin";
	string L=GENOME;L+=".ann.location";
	FILE* BINFILE=File_Open(G.c_str(),"r");
	FILE* Location_File=File_Open(L.c_str(),"r");
	
	fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
	char* Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n"); 
	if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");/*jq: format of genome*/
	char* Marked_Genome=new char[Genome_Size+1];if(!Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 
	
	struct Offset_Record
	{
		char Genome[40];
		unsigned Offset;
	} Temp_OR; 
	
	while (fgets(Temp_OR.Genome,39,Location_File)!=0)//count genomes..
	{
		fgets(Temp_OR.Genome,39,Location_File);
		Genome_CountX++;	
	}
	rewind(Location_File);
	//Offset_Record Genome_Offsets[Genome_CountX+1];

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
		Genome_Count++;	
	}
	Genome_Count--;
	
	char* Split_Point=Org_Genome;//split and write...
	
	//Gene_Hash Genome_List[Genome_Count];
	for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
	{
		String_Hash[Genome_Offsets[i].Genome]=i;
		Genome_List[i].Genome=(Split_Point+=Genome_Offsets[i].Offset);
		Genome_List[i].Index=i;
	}
	
//------------ Analyse read file -------------
	F.Input_File=File_Open(INPUTFILE,"r");
	F.Buffer=(char*) malloc(60000);

	Analyze_File(F);//,L);//Get File type, read length etc..
	printf("Loading C->T Index..");Load_FM_and_Sample(fwfmiCT,revfmiCT,mmPool,CT);printf("[Done]\n");//load FM index structures..
	printf("Loading G->A Index..");Load_FM_and_Sample(fwfmiGA,revfmiGA,mmPool,GA);printf("[Done]\n");
	
	Load_Location(GA.LOCATIONFILE,Annotations);
	Init(F.SOLID,0);//Set SOLiD /Solexa algo for batman..
	LOOKUPSIZE=Get_Lookup_Size(MAX_MISMATCHES,70);//L.STRINGLENGTH);//initialize lookup tables for batman..;
	MLookCT.Lookupsize=MLookGA.Lookupsize=LOOKUPSIZE;
	Build_Tables(fwfmiCT,revfmiCT,MLookCT);
	Build_Tables(fwfmiGA,revfmiGA,MLookGA);

//}-------------------------------- FM File Stuff ---------------------------------------------------
	time_t Start_Time,End_Time;
	time(&Start_Time);

	int Total_Tags=0;

	//Do_Meth(); 
        Thread_Arg T;
	if (THREAD)
	{
		//int Dir_Status=mkdir("batmeth_out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		Launch_Threads(THREAD, Do_Meth,T);
	}
	else
	{
		Do_Meth(NULL);
	}
	Done_Progress();
	time(&End_Time);printf("Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));//Total Reads - %u\n Total_Tags,

}

void *Do_Meth( void *T)
{
	//printf("nihao\n");
	unsigned Actual_Tag=0;
	unsigned Mapped=0;
	READ R,R_CT,R_GA;//+ strand reads..
	READ R_Rev,R_Rev_CT,R_Rev_GA;//Converted reverse strands..
	BATREAD B;
	GUESS G;
	MEMX MF,MC,MF_Rev,MF_GA,MC_GA;
	OUTPUT O_CT,O_GA,O_RevCT,O_RevGA;

	LEN L; //because the range too large,so threads has war.
	L.IGNOREHEAD=0;
	Initial_Length(L);
	if(TRIM_LENGTH) {
		L.STRINGLENGTH_ORG=TRIM_LENGTH;
		printf("Reads are trimmed to %d bp.\n\n",TRIM_LENGTH);
	}
	Split_Read(L.STRINGLENGTH_ORG,L);
	/*LOOKUPSIZE=Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);
	MLookCT.Lookupsize=MLookGA.Lookupsize=LOOKUPSIZE;
        Build_Tables(fwfmiCT,revfmiCT,MLookCT);
        Build_Tables(fwfmiGA,revfmiGA,MLookGA);	
*/
	std::string tmpList1[1000];
	std::string tmpList2[1000];
	std::string tmpList3[1000];
	std::string tmpList4[1000];

	char* Buffer_CT=(char*) malloc(FILEBUFSIZE);
	char* Buffer_RevCT=(char*) malloc(FILEBUFSIZE);
	char* Buffer_GA=(char*) malloc(FILEBUFSIZE);
	char* Buffer_RevGA=(char*) malloc(FILEBUFSIZE);

	O_CT.SAM=0;
	O_CT.PLUSSTRAND=0;
	O_CT.MaxHits=MAXHITS;
	O_GA.MaxHits=MAXHITS;
	O_CT.Offset=0;
	//O_CT.Location_Array=Location_ArrayCT;
	//O_CT.Genome_Offsets=Genome_OffsetsCT;
	O_CT.Genome_Count=Genome_CountCT;
	O_CT.FILETYPE=FQ;
	O_CT.Length_Array[0]=O_CT.Length_Array[1]=O_CT.Length_Array[2]=L.STRINGLENGTH;
	O_RevCT=O_CT;
	O_GA=O_CT;O_RevGA=O_CT;
	O_CT.Buffer=Buffer_CT;//Buffers to write batman output ...
	O_GA.Buffer=Buffer_GA;
	O_RevCT.Buffer=Buffer_RevCT;
	O_RevGA.Buffer=Buffer_RevGA;

	Thread_Arg *TA;
	TA=(Thread_Arg*) T;
	int Thread_ID=TA->ThreadID;

	char ONEFMINDEX=FALSE;

//----------- initialise memory structures ------------------
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	MF.Lookupsize=MF_Rev.Lookupsize=LOOKUPSIZE;
	MF_GA.Lookupsize=MC_GA.Lookupsize=LOOKUPSIZE;
	Copy_MEM(MLookCT,MF,MC,MAX_MISMATCHES);
	Copy_MEM(MLookCT,MF_Rev,MC,MAX_MISMATCHES);
	Copy_MEM(MLookGA,MF_GA,MC_GA,MAX_MISMATCHES);

//----------- Mapping Loop  ------------------
	FILE* Output_File;
	char Out_Name[100];
	if(THREAD==1) {Output_File=File_Open(OUTPUTFILE,"w");Thread_ID=1;}//File_Open(OUTPUTFILE,"a");
	else
	{
		std::string Temp=OUTPUTFILE;
		char Temp_Char[20];sprintf(Temp_Char,".%d",Thread_ID);
		Temp+=Temp_Char;
		Output_File=File_Open(Temp.c_str(),"w");
	}
	//out fq
	std::string tempfq=OUTPUTFILE;tempfq+=".temp.fq";
	TEMP_FQ=File_Open(tempfq.c_str(),"w");
	
	int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
	Init_Progress();
	
	//omp_set_num_threads(THREAD); //Set #threads to use from user
	//#pragma omp parallel default(shared) 
	//{
	while (Read_Tag(F.Input_File,F.FILETYPE,R))
	{	
		if(strlen(R.Tag_Copy)-1!=L.STRINGLENGTH)
		{
			GET_LEN(L,R);
			if(TRIM_LENGTH) {
				L.STRINGLENGTH_ORG=TRIM_LENGTH;
				printf("Reads are trimmed to %d bp.\n\n",TRIM_LENGTH);
			}
			Split_Read(L.STRINGLENGTH_ORG,L);//calculate read portions for Batman algo..
			//L.STRINGLENGTH_ORG=L.STRINGLENGTH;
			//printf("%d\n",L.STRINGLENGTH);
			O_CT.Length_Array[0]=O_CT.Length_Array[1]=O_CT.Length_Array[2]=L.STRINGLENGTH;
        		O_RevCT=O_CT;
        		O_GA=O_CT;O_RevGA=O_CT;
			B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
		}

		if(++Actual_Tag == MAX_READS) break;
		Progress++;
		if (Progress>Number_of_Tags && Thread_ID==1) 
		{
			off64_t Current_Pos=ftello64(F.Input_File);
			off64_t Average_Length=Current_Pos/Actual_Tag+1;//+1 avoids divide by zero..
			Number_of_Tags=(F.File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/F.File_Size);
		}
		int Hits;
		bool Header_Printed = false;
		R_CT=R_Rev_CT=R_Rev_GA=R_GA=R;
		int i;
		//#pragma omp parallel for
		for (i=0;i<=B.StringLength-1;i++){R_Rev.Tag_Copy[B.StringLength-1-i]=R_Rev_GA.Tag_Copy[B.StringLength-1-i]=R_Rev_CT.Tag_Copy[B.StringLength-1-i]=Char_To_C[R.Tag_Copy[i]];};//Reverse complement reads 
		
		int cntHit1 = 0, cntHit2 = 0, cntHit3 = 0, cntHit4 = 0;
		int firstmismatch1 = 0, firstmismatch2 = 0, firstmismatch3 = 0, firstmismatch4 = 0;
		int lastMismatch1 = -1, lastMismatch2 = -1, lastMismatch3 = -1, lastMismatch4 = -1;

		bool Multipe=false;
		minMis=200;
		
		int minTrueMis=100;
		int maxTrueMis=0;
		int SecondMisN=0;
//------------- C->T scan Start ----------------------------------------------
		Hits=0;
		ReplaceCtoT(R_CT.Tag_Copy);
		if(Process_Read(R_CT,B,MF,MC)) //jqjq
		{
			if(Hits=Map_Strand(MAX_MISMATCHES,maxhits,L,fwfmiCT,revfmiCT,MF)){Mapped++;}
			string readString=R.Tag_Copy;
			Print_Hits(SecondMisN,Multipe,minTrueMis,maxTrueMis,'1',Genome_Offsets,readString,MF,L,Output_File,ONEFMINDEX,revfmiCT,fwfmiCT,O_CT);
			*O_CT.Buffer_End='&';*(O_CT.Buffer_End+1)='\n';

			char Banner[2000];
			int D;for(D=0;R.Description[D]!='\n' && R.Description[D]!='\r';D++);R.Description[D]=0;//Make read to C string
			R.Tag_Copy[L.STRINGLENGTH]=0;R.Quality[L.STRINGLENGTH]=0;
			sprintf(Banner,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\n",R.Description,R.Tag_Copy,R.Quality,MF.Stats[0],MF.Stats[1],MF.Stats[2],MF.Stats[3],MF.Stats[4],MF.Stats[5],MF.Stats[6]);
			Get_Records(O_CT.Buffer, tmpList1,firstmismatch1,lastMismatch1,cntHit1,'1');
			

			//------------- G->A reversed scan Start ----------------------------------------------
			Hits=0;
			ReplaceGtoA(R_Rev_GA.Tag_Copy);
			Process_Read(R_Rev_GA,B,MF_GA,MC);
		
			if(Hits=Map_Strand(MAX_MISMATCHES,maxhits,L,fwfmiGA,revfmiGA,MF_GA)){Mapped++;}
			readString=R_Rev.Tag_Copy;
			Print_Hits(SecondMisN,Multipe,minTrueMis,maxTrueMis,'4',Genome_Offsets,readString,MF_GA,L,Output_File,ONEFMINDEX,revfmiGA,fwfmiGA,O_RevGA);
			*O_RevGA.Buffer_End='&';*(O_RevGA.Buffer_End+1)='\n';
			Get_Records(O_RevGA.Buffer, tmpList4,firstmismatch4,lastMismatch4,cntHit4,'4');
			
		
//------------- G->A scan Start ----------------------------------------------
			if(ORIENTATION) 
			{
				Hits=0;
				ReplaceGtoA(R_GA.Tag_Copy);
				Process_Read(R_GA,B,MF_GA,MC);
		
				if(Hits=Map_Strand(MAX_MISMATCHES,maxhits,L,fwfmiGA,revfmiGA,MF_GA)){Mapped++;}
				readString=R.Tag_Copy;
				Print_Hits(SecondMisN,Multipe,minTrueMis,maxTrueMis,'2',Genome_Offsets,readString,MF_GA,L,Output_File,ONEFMINDEX,revfmiGA,fwfmiGA,O_GA);
				*O_GA.Buffer_End='&';*(O_GA.Buffer_End+1)='\n';
				Get_Records(O_GA.Buffer, tmpList2,firstmismatch2,lastMismatch2,cntHit2,'2');
		
//------------- C->T reversed scan Start ----------------------------------------------
				Hits=0;
				ReplaceCtoT(R_Rev_CT.Tag_Copy);
				Process_Read(R_Rev_CT,B,MF,MC);
			
				if(Hits=Map_Strand(MAX_MISMATCHES,maxhits,L,fwfmiCT,revfmiCT,MF)){Mapped++;}
				readString=R_Rev.Tag_Copy;
				Print_Hits(SecondMisN,Multipe,minTrueMis,maxTrueMis,'3',Genome_Offsets,readString,MF,L,Output_File,ONEFMINDEX,revfmiCT,fwfmiCT,O_RevCT);
				*O_RevCT.Buffer_End='&';*(O_RevCT.Buffer_End+1)='\n';
				Get_Records(O_RevCT.Buffer, tmpList3,firstmismatch3,lastMismatch3,cntHit3,'3');

			}
			if(SecondMisN < 300)
			{
				Print_List(Multipe,minTrueMis,Header_Printed,cntHit1,firstmismatch1,lastMismatch1,tmpList1,Banner,Output_File);
				Print_List(Multipe,minTrueMis,Header_Printed,cntHit4,firstmismatch4,lastMismatch4,tmpList4,Banner,Output_File);
				if(ORIENTATION)
				{
					Print_List(Multipe,minTrueMis,Header_Printed,cntHit2,firstmismatch2,lastMismatch2,tmpList2,Banner,Output_File);
					Print_List(Multipe,minTrueMis,Header_Printed,cntHit3,firstmismatch3,lastMismatch3,tmpList3,Banner,Output_File);
				}
			
				if(!Header_Printed)
				{
					flockfile(TEMP_FQ);
					fprintf(TEMP_FQ,"%s\n%s\n+\n%s\n",R.Description,R.Tag_Copy,R.Quality);
					funlockfile(TEMP_FQ);
				}
			}else
			{
				fprintf(Output_File,"@\n%s",Banner);
			}
		}
		
//------------- G->A reversed scan Finish ----------------------------------------------
	}
	return (void *)(Actual_Tag);
}

/* Given a base name, will find file names of related FM Index */
void Build_Names(const char* Genome_Name,FMFILES & FM)//build FM index file names...
{
	int Last_Dash=0;
	char* Name=(char*)Genome_Name;
	for(int i=0;Name[0]!=0;Name++,i++)
	{
		if (Name[0]=='/') 
		{
			Last_Dash=i;
		}
	}
	Name=(char*)Genome_Name+Last_Dash;if(Last_Dash) Last_Dash++;

	char* Command_Line_Buffer=(char*)malloc(5000);
	FM.REVBWTINDEX = (char*)Command_Line_Buffer;
	//if(Last_Dash) Last_Dash=Genome_Name-N+1; else Genome_Name--;
	if(!Last_Dash) Name--;
	strncpy(FM.REVBWTINDEX,Genome_Name,Last_Dash);
	FM.REVBWTINDEX[Last_Dash+0]='r';FM.REVBWTINDEX[Last_Dash+1]='e';FM.REVBWTINDEX[Last_Dash+2]='v';
	strcpy(FM.REVBWTINDEX+Last_Dash+3, Name+1);
	strcat(FM.REVBWTINDEX+Last_Dash+3,".bwt"); 

	FM.BWTFILE=FM.REVBWTINDEX+500;
	strncpy(FM.BWTFILE,Genome_Name,Last_Dash);
	strcpy(FM.BWTFILE+Last_Dash,Name+1);
	strcat(FM.BWTFILE+Last_Dash,".bwt"); 


	FM.REVOCCFILE = FM.BWTFILE+500;
	strncpy(FM.REVOCCFILE,Genome_Name,Last_Dash);
	FM.REVOCCFILE[Last_Dash+0]='r';FM.REVOCCFILE[Last_Dash+1]='e';FM.REVOCCFILE[Last_Dash+2]='v';
	strcpy(FM.REVOCCFILE+Last_Dash+3,Name+1);
	strcat(FM.REVOCCFILE+Last_Dash+3,".fmv"); 


	FM.OCCFILE=FM.REVOCCFILE+500;			
	strncpy(FM.OCCFILE,Genome_Name,Last_Dash);
	strcpy(FM.OCCFILE+Last_Dash,Name+1);
	strcat(FM.OCCFILE+Last_Dash,".fmv"); 

	FM.SAFILE=FM.OCCFILE+500;			
	strncpy(FM.SAFILE,Genome_Name,Last_Dash);
	strcpy(FM.SAFILE+Last_Dash,Name+1);
	strcat(FM.SAFILE+Last_Dash,".sa");

	FM.REVSAFILE = FM.SAFILE+500;
	strncpy(FM.REVSAFILE,Genome_Name,Last_Dash);
	FM.REVSAFILE[Last_Dash+0]='r';FM.REVSAFILE[Last_Dash+1]='e';FM.REVSAFILE[Last_Dash+2]='v';
	strcpy(FM.REVSAFILE+Last_Dash+3,Name+1);
	strcat(FM.REVSAFILE+Last_Dash+3,".sa"); 

	FM.BINFILE=FM.REVSAFILE+500;			
	strncpy(FM.BINFILE,Genome_Name,Last_Dash);
	strcpy(FM.BINFILE+Last_Dash,Name+1);
	strcat(FM.BINFILE+Last_Dash,".bin");

	FM.LOCATIONFILE=FM.BINFILE+500;			
	strncpy(FM.LOCATIONFILE,Genome_Name,Last_Dash);
	strcpy(FM.LOCATIONFILE+Last_Dash,Name+1);
	strcat(FM.LOCATIONFILE+Last_Dash,".ann.location");

	FM.PACFILE=FM.LOCATIONFILE+500;			
	strncpy(FM.PACFILE,Genome_Name,Last_Dash);
	strcpy(FM.PACFILE+Last_Dash,Name+1);
	strcat(FM.PACFILE+Last_Dash,".pac");
}

inline int getMismatch(char* s2) 
{
	int tab=0;
	for (int i=0;s2[i]!='\n';i++)
	{
		if (s2[i]=='\t') { 
			tab++;
			//continue;
		}
		if (tab==4)
		{
			return atoi(s2+i+1);
		}
		else if (tab>4) break;
	}
	return 0;
}
inline int getNtab(char* s2) 
{
	int tab=0;
	for (int i=0;s2[i]!='\n';i++)
	{
		if (s2[i]=='\t') { 
			tab++;
			//continue;
		}
	}
	return tab;
}
bool Get_Records(char* File_Buffer, std::string*  tmpList,int & firstmismatch,int & lastMismatch,int & cntHit,char Count)
{
	int t=0,Skip=0,tab=0;
	char Buffer[MAX_REC_SIZE];
	while(sscanf(File_Buffer,"%[^\n]%n",Buffer,&Skip)) 
	{
		if(Buffer[0] == '&') break;
		tab=getNtab(Buffer);
/*		if(tab==7) 
		{
			File_Buffer+=(Skip+1);
			continue; 
		}
*/		Buffer[0]=Count;//Buffer[Skip]='\n';
		tmpList[t++]=Buffer;
		if (cntHit==0) firstmismatch = getMismatch(Buffer);
		cntHit++;
		if (cntHit==maxhits)lastMismatch = getMismatch(Buffer);
		File_Buffer+=(Skip+1);
	}
	return true;
}

inline void Print_List(bool & Multipe,int & minTrueMis,bool & Header_Printed,int cntHit,int firstmismatch,int lastMismatch,std::string* tmpList,char* Header,FILE *Output_File)
{
	if(MAX_MISMATCHES>3 && Multipe && minTrueMis>3) return;
	flockfile(Output_File);
	//if(cntHit>MAXHITS) cntHit=MAXHITS-1;
	if ((cntHit<maxhits || firstmismatch!=lastMismatch) && cntHit>0) 
	{
		if (!Header_Printed) {fprintf(Output_File,"@\n%s",Header);Header_Printed = true;}
		//if(allHits>=muchHits) cntHit=0;
		int currCnt=0;
		while(cntHit>0) 
		{
			fprintf(Output_File,"%s\n",tmpList[--cntHit].c_str());
			if(currCnt++==MAXHITS) break;
		}
	}
	funlockfile(Output_File);
}

inline void ReplaceCtoT(char* Read)
{
	int i;
	//#pragma omp parallel for
        for (i=0;Read[i]!='\n';i++)
        {
                if (Read[i] == 'C' || Read[i] == 'c') Read[i]='T';
        }
}

inline void ReplaceGtoA(char* Read)
{
	int i;
	//#pragma omp parallel for
        for (i=0;Read[i]!='\n';i++)
        {
                if (Read[i] == 'G' || Read[i] == 'g') Read[i]='A';
        }
}

option Long_Options_Decode[]=
{
	{"help",0,NULL,'h'},
	{"inputfile",1,NULL,'i'},
	{"outputfile",1,NULL,'o'},
	{"genome",1,NULL,'g'},
	{"maxhits",1,NULL,'m'},
	{"maxtags",1,NULL,'t'},
	{"orientation",1,NULL,'O'},
	{"thread",1,NULL,'p'},
	{"forcelength",1,NULL,'F'},
	{0,0,0,0}
};


void Parse_Command_line(int argc, char* argv[],char* & GENOME,char* & CT_GENOME,char* & GA_GENOME,char* & INPUTFILE,char* & OUTPUTFILE,int & MAXHITS, int & MAXMISMATCH,unsigned & MAX_READS)
{
	int Current_Option=0;
	const char* Short_Options ="hi:g:o:n:m:t:O:p:F:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	const char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --inputfile | -i <filename>\t\t Name of input file\n"
" --genome | -g | G->A genome...\n"
" --outputfile | -o <filename>\t\t Name of output file\n"
" --maxmis | -n <integer>\t\t Number of mismatches\n"
" --thread | -p <integer>\t\t Number of CPU threads\n"
" --non_directional | -O (All 4 possible orientations | Only 2 orientations from Original DNA Strands[default])\n"
" --maxhits | -m <INT> (Size for counting list)\n"
" --forcelength | -F <INT> (5' read length to map)\n"
;
	char* Name;int Last_Dash;char* Genome_Name;
	char* Command_Line_Buffer=(char*)malloc(5000);
	GENOME=0;

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options_Decode, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'i':
				INPUTFILE=optarg;
				break;
			case 'o':
				OUTPUTFILE=optarg;
				break;
			/*case 'C':
				CT_GENOME=optarg;
				break;*/
			case 'g':
				GENOME=optarg;
				break;
			/*case 'G':
				GA_GENOME=optarg;
				break;*/
			case 'n':
				MAXMISMATCH=atol(optarg);
				break;
			case 'm':
				MAXHITS=atol(optarg);
				break;
			case 't':
				MAX_READS=atol(optarg);
				break;
			case 'O':
				ORIENTATION=true;
				break;
			case 'p':
				THREAD=atoi(optarg);
				break;
			case 'F':
				TRIM_LENGTH=atoi(optarg);
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
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
                if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);} else Thread_Num++;
        }
        //printf("%d Threads runnning ...\n",Thread_Num);
        pthread_attr_destroy(&Attrib);

        for (int i=0;i<NTHREAD;i++)
        {
                pthread_join(Thread_Info[i].Thread,NULL);
        }
}


