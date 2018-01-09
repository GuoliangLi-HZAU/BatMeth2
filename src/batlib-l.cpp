#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <signal.h> 
#include <assert.h> 
//#include <dvec.h>
#include <getopt.h>
#include "zlib.h"
#include <queue>
#include <map>
#include <ctype.h>

#include "batlib-l.h"
extern int LEAST_MIS;
extern unsigned Genome_Count;
extern int NCOUNT;
extern int GOOD_N_LEN;
extern std::map <unsigned, Ann_Info> Annotations;
char OUTPUT_ZIPPED;
extern char SOLID;
extern char RECOVER_N;
extern char NPOLICY;
extern int TRIM_LENGTH;
extern int Max_MM;
char Random_Array[]="tatacgataggacaatgtcttcgaagcccacgcggtaagccggtcattgcggttgtgcgaacactatcagcctcgctgcatggttaccctgggtggataggacgtttgcccgacattttgacacgcataaaaggtctgtagtgggggtggcacaccataaaccctggggcggctccacgatcgtaaaatcctgcgatctg";
extern const char* Code_To_Char;extern const char* Code_To_CharCAPS;
extern char Char_To_CodeC[];
extern char Char_To_Code[];
extern char Char_To_CharC[];
extern BWT *revfmi,*fwfmi;
extern unsigned Conversion_Factor;
extern char* LOG_SUCCESS_FILE;
extern long File_size;
extern FILE* Log_SFile;
extern int SW_STRING_BUFFER;
extern std::map <std::string,int> String_Hash;
extern Gene_Hash Genome_List[80];
extern int MAX_MISMATCHES;
extern int minMis;
extern int maxhits;
extern int allHits;
extern int indelMis;
extern unsigned char* Original_Text_Ori;
int BATBUF=2000;
//{-----------------------------------  INIT ---------------------------------------------------------
void GET_LEN(LEN & L,READ & R)
{
	//L.STRINGLENGTH=strlen(Read.Tag_Copy);
	L.STRINGLENGTH=0;
	for(;R.Tag_Copy[L.STRINGLENGTH]!=0 && R.Tag_Copy[L.STRINGLENGTH]!='\n' && R.Tag_Copy[L.STRINGLENGTH]!='\r';L.STRINGLENGTH++);
	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Analyze_File
 *  Description:  Analyse a file and try to determine its type..
 * =====================================================================================
 */
//void Analyze_File(char & FILETYPE,char & SOLID,FILE* Input_File,int & TAG_COPY_LEN,char & TAB_SEPERATED,int &  PAIR_LENGTH_LEFT,int & PAIR_LENGTH_RIGHT,LEN & L)
void Analyze_File(INFILE & I,LEN & L)
{
	char Description[MAXDES];
	char Tag_Copy[MAXTAG];
	char Current_Tag[MAXTAG];
	char Quality[MAXTAG];
	char Plus[MAXTAG];
	unsigned Last=0;

	fseek(I.Input_File, 0L, SEEK_END);
	I.File_Size = ftello64(I.Input_File);
	rewind(I.Input_File);

	for(;;)//ignore comments...
	{
		//gzgets(Input_File,Description,MAXDES);
		if(!fgets(Description,MAXDES,I.Input_File)) {if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Analyze_File(): error reading file...\n");printf("Analyze_File(): error reading file...\n");exit(-1);}
		if (Description[0] != '#') break;
		Last=ftello64(I.Input_File);//mark last comment...
	}
	if(!fgets(Current_Tag,MAXTAG,I.Input_File)){if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Analyze_File(): error reading file...\n");printf("Analyze_File(): error reading file...\n");exit(-1);}
	if(Current_Tag[2]>='0' && Current_Tag[2] <='3') I.SOLID=TRUE;else I.SOLID =FALSE;//csfasta has numbers..
	for(I.TAG_COPY_LEN=0;Current_Tag[I.TAG_COPY_LEN]!='\n' && Current_Tag[I.TAG_COPY_LEN]!='\r' && Current_Tag[I.TAG_COPY_LEN]!=0;I.TAG_COPY_LEN++);// Find the length of tag line.should be the same length is quality is present.. 
	for(L.STRINGLENGTH=0;Current_Tag[L.STRINGLENGTH]!='\n' && Current_Tag[L.STRINGLENGTH]!='\r' && Current_Tag[L.STRINGLENGTH]!=0 && Current_Tag[L.STRINGLENGTH]!=PAIR_END_SEPERATOR;L.STRINGLENGTH++);//scan for a split that indicates tab seperated PET

	if(Current_Tag[L.STRINGLENGTH]==PAIR_END_SEPERATOR) 
	{
		I.TAB_SEPERATED=TRUE;//we have pair ended tags..
		if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Analyse_File(): Tab seperated files not supported yet ..\n");printf("Analyse_File(): Tab seperated files not supported yet ..\n");exit(-1);
		I.PAIR_LENGTH_LEFT=L.STRINGLENGTH;
		for(I.PAIR_LENGTH_RIGHT=0;Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\n' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\r' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!=0;I.PAIR_LENGTH_RIGHT++);
	}
	else I.TAB_SEPERATED=FALSE;
	if(!fgets(Quality,MAXTAG,I.Input_File)){if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Analyze_File(): error reading file...\n");printf("Analyze_File(): error reading file...\n");exit(-1);}//plus
	if (Quality[0]=='>') I.FILETYPE=FA;
	else 
	{
		I.FILETYPE=FQ;
		if(Quality[0] != '+' || Description[0] != '@') {if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Init_Variables: Cannot determine file type ...\n");printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
	}
	fseek(I.Input_File,Last,SEEK_SET);//go top
	//if (!TRIM_LENGTH && L.STRINGLENGTH>75) 
/**/	TRIM_LENGTH=75; 
	if (TRIM_LENGTH) L.STRINGLENGTH=TRIM_LENGTH;
	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
	if (I.SOLID) L.IGNOREHEAD=2;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Lookup_Size
 *  Description:  calulate lookupsize for given parameters...
 * =====================================================================================
 */
int Get_Lookup_Size(char MAX_MISMATCHES,char STRINGLENGTH)
{
	int LOOKUPSIZE = 6;
	if (STRINGLENGTH < 28) LOOKUPSIZE =3;
	if (STRINGLENGTH <=51 && MAX_MISMATCHES >5) LOOKUPSIZE = 3;
	return LOOKUPSIZE;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Init
 *  Description:  Initialize general settings..
 * =====================================================================================
 */

void Init(char Solid,char Npolicy)
{
	SOLID=Solid;
	NPOLICY=Npolicy;
	Char_To_Code['N']=0;Char_To_Code['n']=0;Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_Code['0']=0;Char_To_Code['1']=1;Char_To_Code['2']=2;Char_To_Code['3']=3;//for SOLiD
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Allocate_Memory
 *  Description:  Alocate memory for a scan...
 * =====================================================================================
 */
void Allocate_Memory(MEMX & M)
{
	int STRINGLENGTH =36;
	int MAXSTRINGLEN=255;
	int Max_Allocate=1;
	int Max_Limit=5;//3;//5;//3
	for (int i=0;i<Max_Limit-1;i++) Max_Allocate=Max_Allocate*STRINGLENGTH;
	M.BMHStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.FSHStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.FSHStackX0X=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.FSSStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.FSSStackX=(SARange*)malloc(sizeof(SARange)*42*MAXSTRINGLEN);	
	M.BMStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.BMStackX=(SARange*)malloc(sizeof(SARange)*42*MAXSTRINGLEN);	
	M.BMStack_X11=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.BMStack_X11H=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.PSBStack=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	

	M.Hit_Array=(SARange*)malloc(sizeof(SARange)*1000);

	M.Exact_Match_Forward=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);	
	M.Exact_Match_Backward=(SARange*)malloc(sizeof(SARange)*2*MAXSTRINGLEN);
	M.ARRAY_BOUND =sizeof(SARange)*4*Max_Allocate;
	M.Left_Mishits=(SARange*)malloc(M.ARRAY_BOUND);	
	M.Right_Mishits=(SARange*)malloc(M.ARRAY_BOUND);	
	M.Mismatches_Backward=(SARange*)malloc(M.ARRAY_BOUND);//STRINGLENGTH*STRINGLENGTH);//*STRINGLENGTH);	
	M.Mismatches_Forward=(SARange*)malloc(M.ARRAY_BOUND);//STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);	
        M.END_BOUND=sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH;
	M.Two_Mismatches_At_End_Forward=(SARange*)malloc(M.END_BOUND);	
	M.Two_Mismatches_At_End=(SARange*)malloc(M.END_BOUND);	
	M.Possible_20=(SARange*)malloc(M.END_BOUND);///no need extra stringlength.	
	M.Possible_02=(SARange*)malloc(sizeof(SARange)*2*STRINGLENGTH*STRINGLENGTH*STRINGLENGTH);///no need..	
	M.END_BOUND=M.END_BOUND/sizeof(SARange);
	M.ARRAY_BOUND=M.ARRAY_BOUND/sizeof(SARange);

	if (NULL== M.Possible_02||NULL==M.Mismatches_Backward||NULL==M.Two_Mismatches_At_End_Forward||NULL==M.Two_Mismatches_At_End||NULL==M.BMHStack||NULL==M.FSHStackX0X||NULL==M.FSHStack||NULL==M.FSSStack||NULL==M.FSSStackX||NULL==M.BMStack_X11H||NULL==M.BMStack_X11||NULL==M.BMStack||NULL==M.BMStackX||NULL==M.Exact_Match_Backward||NULL==M.Exact_Match_Forward||NULL==M.Exact_Match_Forward||NULL==M.Mismatches_Forward||NULL==M.PSBStack)
	{
		if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Allocate_Memory(): out of memory");
		printf("Allocate_Memory(): out of memory");
		exit(1);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Split_Read
 *  Description:  Calculate the effective string length and Counts the length of halves 
 *  		  of read used be the algorithm..
 * =====================================================================================
 */
void Split_Read(const int STRINGLENGTH, LEN & L)
{

	L.STRINGLENGTH=STRINGLENGTH-L.IGNOREHEAD;
	L.LH=L.STRINGLENGTH/2;//calculate tag portions...
	L.LHQL=L.LH/2;
	if ((L.STRINGLENGTH % 2)) {L.LH++;L.LHQL++;}	
	L.LHQR=L.LH-L.LHQL;
	L.RH=L.STRINGLENGTH-L.LH;
	L.RHQL=L.RH/2;L.RHQR=L.RH-L.RHQL;
//For 10 mismatch extension...
	L.STRINGLENGTHl=L.RH;//calculate tag portions...
	L.LHl=L.STRINGLENGTHl/2;//calculate tag portions...
	L.LHQLl=L.LHl/2;
	if ((L.STRINGLENGTHl % 2)) {L.LHl++;L.LHQLl++;}	
	L.LHQRl=L.LHl-L.LHQLl;
	L.RHl=L.STRINGLENGTHl-L.LHl;
	L.RHQLl=L.RHl/2;L.RHQRl=L.RHl-L.RHQLl;

	L.STRINGLENGTHr=L.LH;//calculate tag portions...
	L.LHr=L.STRINGLENGTHr/2;//calculate tag portions...
	L.LHQLr=L.LHr/2;
	if ((L.STRINGLENGTHr % 2)) {L.LHr++;L.LHQLr++;}	
	L.LHQRr=L.LHr-L.LHQLr;
	L.RHr=L.STRINGLENGTHr-L.LHr;
	L.RHQLr=L.RHr/2;L.RHQRr=L.RHr-L.RHQLr;

}
/*
   Low_QualityC=(char*)malloc(STRINGLENGTH+1);
   Low_QualityF=(char*)malloc(STRINGLENGTH+1);
   Low_Quality=Low_QualityF;

   NLocations=(char*)malloc(STRINGLENGTH+1);

   Do_All=(char*)malloc(STRINGLENGTH+1);
*/

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Copy_MEM 
 *  Description:  Copy Forward and Complement memory structures for scan
 * =====================================================================================
 */
void Copy_MEM(MEMLOOK & M,MEMX & F,MEMX & C,char MAX_MISMATCHES)
{
	/*M.Forward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Forward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));*/

	F.Forward_Start_LookupX=C.Forward_Start_LookupX=M.Forward_Start_LookupX;
	F.Forward_End_LookupX=C.Forward_End_LookupX=M.Forward_End_LookupX;
	F.Backward_Start_LookupX=C.Backward_Start_LookupX=M.Backward_Start_LookupX;
	F.Backward_End_LookupX=C.Backward_End_LookupX=M.Backward_End_LookupX;

	if (MAX_MISMATCHES >5) F.Stat_Size=C.Stat_Size=7; else F.Stat_Size=C.Stat_Size=MAX_MISMATCHES+1;
	Allocate_Memory(F);
	Allocate_Memory(C);

/*	unsigned Write_Buf_Size=((sizeof(Mismatches_Record)+sizeof(Output_Record)+3+100)*1500);//(MAXHITS+1));
	if(!(M.Write_Buffer=(char*)malloc(Write_Buf_Size))) 
	{
		printf("Copy_MEM():out of memory\n");exit(1);
	}
	F.Write_Buffer=C.Write_Buffer=M.Write_Buffer;
*/

}

/*void Copy_MEM(MEM & M,MEMX & F,char MAX_MISMATCHES)
{
	M.Forward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Forward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));

	F.Forward_Start_LookupX=M.Forward_Start_LookupX;
	F.Forward_End_LookupX=M.Forward_End_LookupX;
	F.Backward_Start_LookupX=M.Backward_Start_LookupX;
	F.Backward_End_LookupX=M.Backward_End_LookupX;

	if (MAX_MISMATCHES >5) F.Stat_Size=7; else F.Stat_Size=MAX_MISMATCHES+1;
	Allocate_Memory(F);

	unsigned Write_Buf_Size=((sizeof(Mismatches_Record)+sizeof(Output_Record)+3+100)*1500);//(MAXHITS+1));
	if(!(M.Write_Buffer=(char*)malloc(Write_Buf_Size))) 
	{
		printf("Copy_MEM():out of memory\n");exit(1);
	}

}*/
//}-----------------------------------  INIT ---------------------------------------------------------

//{----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Indexes
 *  Description:  Load reverse and forward indexes...
 * =====================================================================================
 */


void Load_Indexes(BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,FMFILES & FM)//char *BWT,char *OCC,char *SA,char *REVBWT,char *REVOCC,char *REVSA)
{
	fwfmi=initFMI(FM.BWTFILE,FM.OCCFILE,FM.SAFILE,mmPool);//Load FM indexes
	revfmi=initFMI(FM.REVBWTINDEX,FM.REVOCCFILE,FM.REVSAFILE,mmPool);
	unsigned SOURCELENGTH = fwfmi->textLength;
	if (SOURCELENGTH!=revfmi->textLength)
	{ 
		if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Load_Indexes():FM index load error \n"); 
		printf("Load_Indexes():FM index load error \n"); 
		exit(-1);
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initFMI
 *  Description:  Opens FM index fmiFile
 * =====================================================================================
 */
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName,const char* SAFile, MMPool *mmPool) 
{
	BWT *fmi;
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	//fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL);//Load FM index
	fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SAFile, NULL, NULL, NULL);//Load FM index
	return fmi;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_SARange
 *  Description:  gets the SA range of strings having prefix [New_Char][Range]
 * =====================================================================================
 */

SARange Get_SARange( char New_Char,struct SARange Range,BWT *fmi)
{

	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		Range.Start=0;
	}
	return Range;

}

void Get_SARange_Fast( char New_Char,struct SARange & Range,BWT *fmi)
{

	Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.Start, New_Char) + 1;
	Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Range.End+1, New_Char);
	if (Range.End<Range.Start) 
	{
		//Range.End=0;
		Range.Start=0;
	}
}

void Get_SARange_Fast_2( long New_Char,struct SARange & Start_Range, struct SARange & Dest_Range,BWT *fmi)
{

	Dest_Range.Start = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Start_Range.Start, New_Char) + 1;
	Dest_Range.End = fmi->cumulativeFreq[New_Char] + BWTOccValue(fmi, Start_Range.End+1, New_Char);
	if (Dest_Range.End<Dest_Range.Start) 
	{
		//Range.End=0;
		Dest_Range.Start=0;
	}
}
//}----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------


//{----------------------------------- FILE HANDLING ---------------------------------------------------------

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
	if(strcmp(File_Name,"-")==0)
	{
		return stdout;
	}	
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		if (LOG_SUCCESS_FILE) fprintf(Log_SFile,"File %s Cannot be opened ....\n",File_Name);
		printf("File %s Cannot be opened ....\n",File_Name);
		exit(-1);
	}
	else return Handle;
}

/*void File_OpenZ(const char* File_Name,const char* Mode,gzFile & Handle)
{
	Handle=gzopen(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
}*/

size_t fwriteX ( const void * ptr, size_t size, size_t count, void* stream )
{
	/*if (OUTPUT_ZIPPED)
	{
		return gzwrite ( (gzFile*) stream ,ptr, size*count);
	}
	else*/
	{
		return fwrite ( ptr, size, count, (FILE*) stream );
	}
}


void fprintfX(void* Handle,char* Format, char* String)
{

	/*if (OUTPUT_ZIPPED)
	{
		gzprintf((gzFile*) Handle, Format, String);
	}
	else*/
	{
		fprintf((FILE*) Handle, Format, String);
	}
}

//}----------------------------------- FILE HANDLING ---------------------------------------------------------
std::string Read_Sam(FILE *Input_File,const char FILETYPE, READ & Read )
{
	flockfile(Input_File);
	char s2t[BATBUF];char buf2[MAXDES];char buf1[MAXDES];
	char Tag_Copy[MAXDES];char Quality[MAXDES];
	if (fgets(s2t,BATBUF,Input_File)!= 0)// read online fgets(n,MAXDES,Input_File)
	{
		sscanf(s2t,"%s%*d%*s%*d%*d%*s%*s%*d%*d%s%s",Read.Description,Tag_Copy,Quality);
		char s2tSplit[BATBUF];strcpy(s2tSplit,s2t);
		
		    const char *d = "\t";
		    char *p;
		    p = strtok(s2tSplit,d);
		    while(p)
		    {
		    	 if(p[0]=='Y' && p[1]=='C' )
		        	sscanf(p,"YC:Z:%s",buf1);// CT GA
		        else if(p[0]=='Y' && p[1]=='D')
		        	sscanf(p,"YD:Z:%s",buf2);//f r
		        p=strtok(NULL,d);
		    }
		if( (!strcmp(buf2,"r") && !strcmp(buf1,"CT")) || ( !strcmp(buf2,"f") && !strcmp(buf1,"GA") ) )
		{
			int Lens=0;
			for(;Tag_Copy[Lens]!=0 && Tag_Copy[Lens]!='\n' && Tag_Copy[Lens]!='\r';Lens++);
			for (int i=0;i<=Lens-1;i++)
			{
				Read.Tag_Copy[Lens-1-i]=Char_To_CharC[Tag_Copy[i]];
				Read.Quality[Lens-1-i]=Quality[i];
			};
		}else
		{
			strcpy(Read.Tag_Copy,Tag_Copy);
			strcpy(Read.Quality,Quality);
		}
		char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
		funlockfile(Input_File);
		return s2t;
	}
	else 
	{
		funlockfile(Input_File);
		std::string False="FALSE";
		return False;
	}
}

//{-----------------------------------  READ ROUTINES ---------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Tag
 *  Description:  Read a line from FASTA/FASTQ
 *  		  return true if successfull...
 * =====================================================================================
 */
char Read_Tag(FILE *Input_File,const char FILETYPE, READ & Read )
{
	flockfile(Input_File);
	if (fgets(Read.Description,MAXDES,Input_File)!=0)// read a tag...
	{
		char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
		//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
		if(!fgets(Read.Tag_Copy,MAXDES,Input_File)) {if(LOG_SUCCESS_FILE) fprintf (Log_SFile,"Read_Tag():Error reading file..\n"); printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
		if (FILETYPE == FQ)
		{
			//gzgets(Input_File,Plus,MAXTAG);//plus
			if(!fgets(Read.Plus,MAXTAG,Input_File)){if(LOG_SUCCESS_FILE) fprintf (Log_SFile,"Read_Tag():Error reading file..\n"); printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
			//gzgets(Input_File,Quality,MAXTAG);//phred
			if(!fgets(Read.Quality,MAXTAG,Input_File)){if(LOG_SUCCESS_FILE) fprintf (Log_SFile,"Read_Tag():Error reading file..\n"); printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
		}
		else
		{
			Read.Quality[0]='*';Read.Quality[1]=0;
		}
		funlockfile(Input_File);
		return TRUE;
	}
	else 
	{
		funlockfile(Input_File);
		return FALSE;
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Process_Read
 *  Description:  Convert a read to binary and find its rev complement...
 * 		  Load B.StringLength, B contains converted reads.
 * =====================================================================================
 */
bool Process_Read(READ & RawR,READ & R, BATREAD & B,MEMX & MF, MEMX & MC)
{
	int j=0;
	static int Random_Pointer=0;
	int NCumulative[256];

	R.NCount=0;
	for (unsigned i=0;i<=B.StringLength-1;i++)
	{
		if (R.Tag_Copy[i+B.IGNOREHEAD] == 'n' || R.Tag_Copy[i+B.IGNOREHEAD]=='N')
		{
			R.NCount++;
			if (NPOLICY)
			{
				R.N[j++]=i;R.NLocations[i]=TRUE;
				R.Tag_Copy[i+B.IGNOREHEAD]=Random_Array[Random_Pointer++];R.N[j++]=R.Tag_Copy[i+B.IGNOREHEAD];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}else
			{
				R.N[j++]=i;R.NLocations[i]=TRUE;R.N[j++]='N';
			}
		}
		else R.NLocations[i]=FALSE;
		if (RECOVER_N) NCumulative[i]=R.NCount;

		B.Forward[i]=Char_To_Code[R.Tag_Copy[i+B.IGNOREHEAD]];
		B.Forward_raw[i]=Char_To_Code[RawR.Tag_Copy[i+B.IGNOREHEAD]];
		
		if (SOLID) 
		{
			B.Complement[B.StringLength-1-i]=B.Forward[i];
			B.Complement_raw[B.StringLength-1-i]=B.Forward_raw[i];
		}
		else 
		{
			B.Complement[B.StringLength-1-i]=Char_To_CodeC[B.Forward[i]];
			B.Complement_raw[B.StringLength-1-i]=Char_To_CodeC[B.Forward_raw[i]];
		}
	}
//	if(R.NCount>2) {return false;}
	B.Forward[B.StringLength]='+'; B.Forward_raw[B.StringLength]='+';
	B.Complement[B.StringLength]='-'; B.Complement_raw[B.StringLength]='-';
	if (RECOVER_N && R.NCount>NCOUNT && B.StringLength >=75 )//search for best N-less read 
	{
		int S= -1,E,Max=0;
		NCumulative[B.StringLength]=100;//sentinel...
		if(NCumulative[GOOD_N_LEN-1]<1) 
		{
			int j=GOOD_N_LEN-1;
			while(NCumulative[j]<1) {j++;}
			Max=j;S=0;E=j-1;
		}
		for (unsigned i=0;i<=B.StringLength-1-GOOD_N_LEN && B.StringLength-1-GOOD_N_LEN>=0;i++)
		{
			int NC=NCumulative[i+GOOD_N_LEN]-NCumulative[i];
			if (NC<1)
			{
				int j=i+GOOD_N_LEN;
				while(j < B.StringLength && (NCumulative[j]-NCumulative[i])<1 ){j++;}
				if (Max < (j-i-1)) {Max=(j-i-1);S=i+1;E=j-1;}
			}
		}
		if (S>=0) 
		{
			R.N[0]=S;R.N[1]=E;//for(int k=S;k<=E;k++){putchar(R.Tag_Copy[k]);}printf("\n");}
			int j=0;
			for(int i=S;i<=E;i++,j++)
			{
				B.Forward[j]=B.Forward[i];
				B.Complement[j]=B.Complement[B.StringLength-1-E+j];
				
				B.Forward_raw[j]=B.Forward_raw[i];
				B.Complement_raw[j]=B.Complement_raw[B.StringLength-1-E+j];
				//putchar('0'+B.Complement[j]);
				//putchar('0'+B.Forward[j]);
			}
			//printf("\n");
			B.Forward[j]='+';B.Complement[j]='-';
			B.Forward_raw[j]='+';B.Complement_raw[j]='-';
			R.N[2]=j;
#ifndef NDEBUG
			for(int i=0;i<Max;i++)
			{
				assert(B.Complement[Max-1-i]==Char_To_CodeC[B.Forward[i]]);
			}
#endif
		}
		else
		{
			R.N[0]=0;R.N[1]=0;R.N[2]=0;//for(int k=S;k<=E;k++){putchar(R.Tag_Copy[k]);}printf("\n");}
		}
	}


	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Read=R;MF.Current_Tag=B.Forward; MF.Current_Tag_raw=B.Forward_raw;  MF.Hit_Array[0].Start=0;
	MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Read=R;MC.Current_Tag=B.Complement; MC.Current_Tag_raw=B.Complement_raw;  MC.Hit_Array[0].Start=0;
	for (int i=0;i<7;i++) {MF.Stats[i]=0;MC.Stats[i]=0;}
	return true;
}

/*
void Process_Read(READ & R, BATREAD & B,MEMX & MF, MEMX & MC)
{
	int j=0;
	static int Random_Pointer=0;
	int NCumulative[256];

	R.NCount=0;
	for (unsigned i=0;i<=B.StringLength-1;i++)
	{
		if (R.Tag_Copy[i+B.IGNOREHEAD] == 'n' || R.Tag_Copy[i+B.IGNOREHEAD]=='N')
		{
			R.NCount++;
			if (NPOLICY)
			{
				R.N[j++]=i;R.NLocations[i]=TRUE;
				R.Tag_Copy[i+B.IGNOREHEAD]=Random_Array[Random_Pointer++];R.N[j++]=R.Tag_Copy[i+B.IGNOREHEAD];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
		}
		else R.NLocations[i]=FALSE;
		if (RECOVER_N) NCumulative[i]=R.NCount;

		B.Forward[i]=Char_To_Code[R.Tag_Copy[i+B.IGNOREHEAD]];
		if (SOLID) B.Complement[B.StringLength-1-i]=B.Forward[i];
		else B.Complement[B.StringLength-1-i]=Char_To_CodeC[B.Forward[i]];
	}
	B.Forward[B.StringLength]='+';
	B.Complement[B.StringLength]='-';
	if (RECOVER_N && R.NCount>NCOUNT && B.StringLength >=75)//search for best N-less read
	{
		int S= -1,E,Max=0;
		NCumulative[B.StringLength]=100;//sentinel...
		if(NCumulative[GOOD_N_LEN-1]<1) 
		{
			int j=GOOD_N_LEN-1;
			while(NCumulative[j]<1) {j++;}
			Max=j;S=0;E=j-1;
		}
		for (unsigned i=0;i<=B.StringLength-1-GOOD_N_LEN && B.StringLength-1-GOOD_N_LEN>0;i++)
		{
			int NC=NCumulative[i+GOOD_N_LEN]-NCumulative[i];
			if (NC<1)
			{
				int j=i+GOOD_N_LEN;
				while(j < B.StringLength && NCumulative[j]-NCumulative[i]<1 ){j++;}
				if (Max < (j-i-1)) {Max=(j-i-1);S=i+1;E=j-1;}
			}
		}
		if (S>=0) 
		{
			R.N[0]=S;R.N[1]=E;//for(int k=S;k<=E;k++){putchar(R.Tag_Copy[k]);}printf("\n");}
			int j=0;
			for(int i=S;i<=E;i++,j++)
			{
				B.Forward[j]=B.Forward[i];
				B.Complement[j]=B.Complement[B.StringLength-1-E+j];
				//putchar('0'+B.Complement[j]);
				//putchar('0'+B.Forward[j]);
			}
			//printf("\n");
			B.Forward[j]='+';B.Complement[j]='-';
			R.N[2]=j;
#ifndef NDEBUG
			for(int i=0;i<Max;i++)
			{
				assert(B.Complement[Max-1-i]==Char_To_CodeC[B.Forward[i]]);
			}
#endif
		}
		else
		{
			R.N[0]=0;R.N[1]=0;R.N[2]=0;//for(int k=S;k<=E;k++){putchar(R.Tag_Copy[k]);}printf("\n");}
		}
	}


	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Read=R;MF.Current_Tag=B.Forward;MF.Hit_Array[0].Start=0;
	MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Read=R;MC.Current_Tag=B.Complement;MC.Hit_Array[0].Start=0;
	for (int i=0;i<7;i++) {MF.Stats[i]=0;MC.Stats[i]=0;}

}
*/

bool Process_Read_bat(READ & R, BATREAD & B,MEMX & MF, MEMX & MC)
{
	int j=0;
	static int Random_Pointer=0;

	R.NCount=0;
	for (unsigned i=0;i<=B.StringLength-1;i++)
	{
		if (R.Tag_Copy[i+B.IGNOREHEAD] == 'n' || R.Tag_Copy[i+B.IGNOREHEAD]=='N')
		{
			R.NCount++;
			if (NPOLICY)
			{
				R.N[j++]=i;R.NLocations[i]=TRUE;
				R.Tag_Copy[i+B.IGNOREHEAD]=Random_Array[Random_Pointer++];R.N[j++]=R.Tag_Copy[i+B.IGNOREHEAD];
				if (Random_Pointer==sizeof(Random_Array)-1) Random_Pointer=0; 
			}
			else
			{
				R.N[j++]=i;R.NLocations[i]=TRUE;R.N[j++]='N';
			}
		}
		else R.NLocations[i]=FALSE;
		B.Forward[i]=Char_To_Code[R.Tag_Copy[i+B.IGNOREHEAD]];
		if (SOLID) B.Complement[B.StringLength-1-i]=B.Forward[i];
		else B.Complement[B.StringLength-1-i]=Char_To_CodeC[B.Forward[i]];
	}
	if(R.NCount>5) {/*if(PRINT_MISHITS) fprintf(Mishit_File,"%s%s%s%s", Description,TAg_Copy,Plus,Quality);*/ return false;} //jqjq_Hardcode
	//if (NISMISMATCH && MAX_MISMATCHES <= R.NCount ) continue;//Too many N's than mismatches...
	B.Forward[B.StringLength]='+';
	B.Complement[B.StringLength]='-';


	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Read=R;MF.Current_Tag=B.Forward;
	MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Read=R;MC.Current_Tag=B.Complement;
	for (int i=0;i<7;i++) {MF.Stats[i]=0;MC.Stats[i]=0;}
	return true;

}

void Recalibrate_Read(READ & R, BATREAD & B,MEMX & MF, MEMX & MC,LEN & L)
{
	assert(R.NCount);
	B.StringLength=R.N[2];
	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Read=R;MF.Current_Tag=B.Forward;MF.Hit_Array[0].Start=0;
	MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Read=R;MC.Current_Tag=B.Complement;MC.Hit_Array[0].Start=0;
	for (int i=0;i<7;i++) {MF.Stats[i]=0;MC.Stats[i]=0;}

	L.IGNOREHEAD=0;
	Split_Read(B.StringLength,L);

}
//}-----------------------------------  READ ROUTINES ---------------------------------------------------------


//{-----------------------------------  BUILD TABLES ---------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Build_Tables
 *  Description:  build a table to cache first LOOKUPSIZE bases..
 * =====================================================================================
 */
void Build_Tables(BWT *fwfmi, BWT *revfmi,MEMLOOK & M)
{
	Range LRange;
	M.Forward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Forward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_Start_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));
	M.Backward_End_LookupX=(unsigned*)malloc(sizeof(unsigned)*(2<<2*M.Lookupsize));

	for( int i=0;i<4;i++)//fill lookup tables for first character...
	{
		LRange.Start=revfmi->cumulativeFreq[i] + 1;
		LRange.End=revfmi->cumulativeFreq[i + 1];
		LRange.Label=i;
		if(LRange.Start>LRange.End) LRange.Start=0; 
		for(int j=0;j<4;j++) Build_Preindex_Forward(revfmi,LRange, 2, j,M);
		LRange.Start=fwfmi->cumulativeFreq[i] + 1;
		LRange.End=fwfmi->cumulativeFreq[i + 1];
		LRange.Label=i;
		if(LRange.Start>LRange.End) LRange.Start=0; 
		for(int j=0;j<4;j++) Build_Preindex_Backward(fwfmi,LRange, 2, j,M);
	}
}

void Build_Preindex_Backward(BWT* fwfmi,Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if(Range.Start)
		{
			Range.Start= fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
			Range.End= fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		M.Backward_Start_LookupX[Range.Label]=Range.Start;
		M.Backward_End_LookupX[Range.Label]=Range.End;
		assert(M.Backward_Start_LookupX[Range.Label]<=M.Backward_End_LookupX[Range.Label]);
	}
	else
	{


		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		if(Range.Start)
		{
			Range.Start = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
			Range.End = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Level ++;

		for ( int i=0;i<4;i++)
		{
			Build_Preindex_Backward(fwfmi, Range, Level,i,M);
		}

	}

}

void Build_Preindex_Forward(BWT *revfmi, Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if(Range.Start)
		{
			Range.Start= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
			Range.End= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		M.Forward_Start_LookupX[Range.Label]=Range.Start;
		M.Forward_End_LookupX[Range.Label]=Range.End;
		assert(M.Forward_Start_LookupX[Range.Label] <=M.Forward_End_LookupX[Range.Label]);
	}
	else
	{


		assert(Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		if (Range.Start)
		{
			Range.Start = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
			Range.End = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		Level ++;
		for ( int i=0;i<4;i++)
		{
			Build_Preindex_Forward( revfmi,Range, Level,i,M);
		}

	}

}
//}-----------------------------------  BUILD TABLES ---------------------------------------------------------

//{-----------------------------------  SCAN ROUTINES ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Zero_Mismatch
 *  Description:  scans for zero mismatch occurrences ...
 *  		  Hits are from M.Hit_Array[0..M.Hit_Array_Ptr-1];
 *  		  returns hits found...
 * =====================================================================================
 */

unsigned Zero_Mismatch(bool batmeth1,char* Current_Tag,LEN L, BWT* revfmi,MEMX & M,BWT* fwfmi)
{
	int c;
	SARange Range;
	M.Left_Mishits_Pointer=0;
	M.Right_Mishits_Pointer=0;
	M.Possible_20_Pointer=0;
	M.Possible_02_Pointer=0;
	M.Mismatches_Forward_Pointer=0;
	M.Mismatches_Backward_Pointer=0;
	M.Two_Mismatches_At_End_Pointer=0;
	M.Two_Mismatches_At_End_Forward_Pointer=0;
	M.FMIndex=REVERSE;

	if(M.Lookupsize ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Range.Start=M.Forward_Start_LookupX[c];Range.End=M.Forward_End_LookupX[c];
	Range.Level=M.Lookupsize+1;Range.Mismatches=0;Range.Skip=0;
	Range.Mismatch_Char=0;
	Search_Forwards_Exact(batmeth1,Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,M,L.LH,revfmi,fwfmi);//Find exact matches and report... if not found get the range for 0|?
	return M.Hits;

}


//{------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------
unsigned One_Mismatch(bool batmeth1,char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{

	M.FMIndex=REVERSE;
	SARange Range;
	int c;

	Range=M.Exact_Match_Forward[L.LH-1];//[Start+L.LH];
	if(Range.Start )//&& Range.Tag == Actual_Tag)//if there are hits of the form 0|?
	{
		Range.Level=1;
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_Forwards(Range,1,LH+1,RH,revfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
			if(MAXHITS==Hits) continue;
		}*/

		//Do_Branch=Do_All;
		Search_Forwards(batmeth1,Current_Tag,Range,1,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
		if(MAXHITS<=M.Hits) return M.Hits;


	}		
	M.FMIndex=FORWARD;

	if(M.Lookupsize ==3)
	{
		c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4) | (Current_Tag[L.STRINGLENGTH-1-3]<<6) | Current_Tag[L.STRINGLENGTH-1-4]<<8 | (Current_Tag[L.STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	Range.Start=M.Backward_Start_LookupX[c];Range.End=M.Backward_End_LookupX[c];Range.Level=M.Lookupsize+1;
	Range.Mismatches=0;Range.Skip=0;Range.Mismatch_Char=0;
	
	Search_Backwards_Exact( Current_Tag,Range,L.STRINGLENGTH,L.RH,fwfmi,M);//Backward scan for ?|0

	if(Range.Start)//if there are possible hits of the form ?|0
	{
		Range.Level=1;
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			Search_Backwards(Range,1,LH,LH,fwfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
			if(MAXHITS<=Hits) return;
		}*/

		//Do_Branch=Do_All;
		Search_Backwards(batmeth1,Current_Tag,Range,1,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	return M.Hits;

}

void Search_Backwards(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,MEMX & M,BWT *revfmi)
{
	if (!Tag.Start) return;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	int BMStack_Top=0;
	M.BMStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMStack[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,fmi,M,revfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches==Count)
							{
								if(batmeth1) Print_LocationX(Temp_Range,M);
								else Print_LocationX(Temp_Range,M,revfmi,fmi);
								if(MAXHITS<=M.Hits) return;
							}
							else continue;
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMStack[BMStack_Top]=Temp_Range;
						}
					}
					else 
					{
						if(5 > Count)// 2 mismatches...
						{
							if(Temp_Range.Level != StringLength) Temp_Range.Level++;
							//if((Start-Temp_Range.Level) != 0) Temp_Range.Level++;
							else // 2 mismatches with the last at the end?
							{
								if(M.Two_Mismatches_At_End_Pointer < M.END_BOUND)
								{
									M.Two_Mismatches_At_End[M.Two_Mismatches_At_End_Pointer]=Temp_Range;
									M.Two_Mismatches_At_End_Pointer++;
								}
								continue;
							}
							if (M.Mismatches_Backward_Pointer < M.ARRAY_BOUND)
							{
								M.Mismatches_Backward[M.Mismatches_Backward_Pointer]=Temp_Range;
								M.Mismatches_Backward_Pointer++;
							}
						}
						continue;
					}
				} 
			}
		}
	}
	return;
}

void Search_Backwards_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M,BWT *revfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					if(batmeth1) Print_LocationX(Tag,M);
					else Print_LocationX(Tag,M,revfmi,fmi);
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else 
		{
			if(5 >= Tag.Mismatches && 5 > Count)// 2 mismatches
			{
				if(!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else//two mismatches with the last at the end ... 
				{
					//if(Tag.Skip) Tag.Start=Tag.End;
					if(M.Two_Mismatches_At_End_Pointer < M.END_BOUND)
					{
						M.Two_Mismatches_At_End[M.Two_Mismatches_At_End_Pointer]=Tag;
						M.Two_Mismatches_At_End_Pointer++;
					}
					return;
				}
				if(M.Mismatches_Backward_Pointer < M.ARRAY_BOUND)
				{
					M.Mismatches_Backward[M.Mismatches_Backward_Pointer]=Tag;
					M.Mismatches_Backward_Pointer++;
				}
			}
			return;
		} 
	}
}

void Branch_Detect_Backwards (const char* Current_Tag,const struct SARange Tag,BWT *fmi,int Start,unsigned* Branch_Characters,SARange *Branch_Ranges)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		assert(1==2);
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else */if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Start-Tag.Level] && Branch != Current_Tag[Start-Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}
//}------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------

unsigned Two_Mismatch(bool batmeth1,char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	int c;
	SARange Range,TRange;
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			(M.Two_Mismatches_At_End_Forward[i]).Mismatch_Pos[1]= (L.STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M,revfmi,fwfmi);
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Forward_Pointer=0;

	M.FMIndex=FORWARD;
	if(M.Two_Mismatches_At_End_Pointer)
	{
		//Do_Branch=Do_All;
		for(int i=0;i<M.Two_Mismatches_At_End_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End[i],M,revfmi,fwfmi);//Mismatches of the form 2|0, with one mismatch at the first position
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}

	M.Two_Mismatches_At_End_Pointer=0;
	M.FMIndex=REVERSE;
	M.Possible_03_Pointer=M.Mismatches_Forward_Pointer;
	if(M.Mismatches_Forward_Pointer)
	{
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_03_Pointer-1;i>=0;i--)
			{
				Search_Forwards(Mismatches_Forward[i],2,LH+1,RH,revfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) continue;
		}

		Do_Branch=Do_All;*/
		for(int i=M.Possible_03_Pointer-1;i>=0;i--)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],2,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=FORWARD;
	M.Possible_30_Pointer=M.Mismatches_Backward_Pointer;
	if(M.Mismatches_Backward_Pointer)
	{
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_30_Pointer-1;i>=0;i--)
			{
				Search_Backwards(Mismatches_Backward[i],2,LH,LH,fwfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) continue;
		}

		Do_Branch=Do_All;*/
		for(int i=M.Possible_30_Pointer-1;i>=0;i--)
		{
			Search_Backwards(batmeth1,Current_Tag,M.Mismatches_Backward[i],2,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	//----------------------------------------------------------------------------------------------------------------------------------------
	if(M.Lookupsize==3)
	{
		c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4) | (Current_Tag[L.STRINGLENGTH-1-3]<<6) | Current_Tag[L.STRINGLENGTH-1-4]<<8 | (Current_Tag[L.STRINGLENGTH-1-5]<<10);//Use lookup table...
	}
	Range.Start=M.Backward_Start_LookupX[c];Range.End=M.Backward_End_LookupX[c];
	Range.Mismatches=0;Range.Level=M.Lookupsize+1;
	Range.Mismatch_Char=0;Range.Skip=0;

	Search_Backwards_Exact_X0( Current_Tag,Range,L.STRINGLENGTH,L.RHQR,fwfmi);// ?|?|0
	Range.Level=1;

	/*if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		Search_Backwards_X10(Range,1,LH + RHQL, RHQL,fwfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) continue;
	}*/

	//Do_Branch=Do_All;
	Search_Backwards_X10(batmeth1,Current_Tag,Range,1,L.LH + L.RHQL, L.RHQL,MAXHITS,L,fwfmi,M,revfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if(MAXHITS<=M.Hits) return M.Hits;
	//----------------------------------------------------------------------------------------------------------------------------------------
	if(M.Lookupsize==3)
	{
		c=Current_Tag[L.LH+0] | (Current_Tag[L.LH+1]<<2) | (Current_Tag[L.LH+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[L.LH+0] | (Current_Tag[L.LH+1]<<2) | (Current_Tag[L.LH+2]<<4) | (Current_Tag[L.LH+3]<<6) | Current_Tag[L.LH+4]<<8 | (Current_Tag[L.LH+5]<<10);//Use lookup table...
	}
	Range.Start=M.Forward_Start_LookupX[c];Range.End=M.Forward_End_LookupX[c];
	Range.Mismatches=0;Range.Level=M.Lookupsize+1; Range.Skip=0;
	Range.Mismatch_Char=0;
	Search_Forwards_0X(Current_Tag,Range,L.LH+1,L.RHQL,revfmi);
	Range.Level=1;TRange=Range;
	/*if(USEQUALITY)
	{
		Do_Branch=Low_Quality;

		Search_X01(Range,1,LH + RHQL +1,RHQR,revfmi);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
		if(MAXHITS==Hits) continue;
	}

	Do_Branch=Do_All;*/
	Search_X01(batmeth1,Current_Tag,TRange,1,L.LH + L.RHQL +1,L.RHQR,MAXHITS,L,fwfmi,revfmi,M);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	return M.Hits;
}

void Search_X01(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT* fwfmi,BWT *revfmi,MEMX & M)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	M.FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_X01_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,L,fwfmi,revfmi,M);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Reverse(Current_Tag,Temp_Range,L.STRINGLENGTH,L.RH,L,fwfmi,M);
								if(Temp_Range.Start)
								{
									Temp_Range.Level=1;
									//Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
									Search_Backwards(batmeth1,Current_Tag,Temp_Range,2,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);
									if(MAXHITS<=M.Hits) return;
									//Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
								}
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(M.Possible_02_Pointer < M.END_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							M.Possible_02[M.Possible_02_Pointer]=Temp_Range;
							M.Possible_02_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}

void Search_X01_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fwfmi,BWT *revfmi,MEMX & M)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
		Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					Reverse(Current_Tag,Tag,L.STRINGLENGTH,L.RH,L,fwfmi,M);
					if(Tag.Start)
					{
						Tag.Level=1;
						Search_Backwards(batmeth1,Current_Tag,Tag,2,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);
					}
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if(M.Possible_02_Pointer < M.END_BOUND)
			{
				if(!Tag.Skip)Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				M.Possible_02[M.Possible_02_Pointer]=Tag;
				M.Possible_02_Pointer++;
			}
			return;
		}
	}
}

void Search_Forwards_0X(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{

	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offsets...
	for(;;)	
	{
		Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
		if (Tag.Start!=0)
		{
			if(Tag.Level== StringLength)
			{
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			return;
		}
	}
}

void Search_Backwards_X10(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fwfmi,MEMX & M,BWT *revfmi)
{
	if (!Tag.Start) return;
	int BMHStack_Top=0;
	M.BMHStack[0]=Tag;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	struct SARange Range,Temp_Range;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_X10_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,L,fwfmi,M,revfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fwfmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								Temp_Range.Level=1;
								//Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
								//memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
								Search_Backwards(batmeth1,Current_Tag,Temp_Range,2,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);
								if(MAXHITS<=M.Hits) return;
								//Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
								//memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if (M.Possible_20_Pointer < M.END_BOUND)
						{
							M.Possible_20[M.Possible_20_Pointer]=Temp_Range;
							M.Possible_20_Pointer++;
						}
					}
				} 
			}
		}
	}
	return;
}


void Search_Backwards_X10_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fmi,MEMX & M,BWT *revfmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)//only one mismatch allowed here...
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Backwards(batmeth1,Current_Tag,Tag,2,L.LH,L.LH,MAXHITS,fmi,M,revfmi);
					return;
				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else
		{
			if(M.Possible_20_Pointer < M.END_BOUND)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				M.Possible_20[M.Possible_20_Pointer]=Tag;
				M.Possible_20_Pointer++;
			}
			return;
		}
	}
}

void Search_Backwards_Exact_X0(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
////can optimise..	
	if (!Tag.Start) return;
	for(;;)
	{
		Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fmi);
		if (Tag.Start!=0)
		{
			if(Tag.Level== StringLength)
			{
				return;
			}
			else Tag.Level++;
		} 
		else
		{
			return;//No hit
		}
	}
}

void Search_Exact(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{
	int Level;
	unsigned long Index,Now,First,Last;
	unsigned Branch_Characters[4];
	if (!Tag.Start) return;

	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}
			Level=Tag.Level;
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start+Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					if(Level== StringLength)
					{
						if(!Tag.Skip) Tag.End=Tag.Start; 
						return;	
					}
					else {Level++;continue;}
				} 
				else//mismatch...
				{
					Tag.Start=0;//Tag.End=0;
					return;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					Tag.Start=0;
				}
			} 
			else
			{
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
			}

			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else//Mismatch
			{
				return;
			}

		}
	}
}

void Backwards(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *revfmi,MEMX & M)
{	

	unsigned Temp=0;
	char Mismatch_Count=Tag.Mismatches;
	unsigned Branch_Characters[4];
	int Temp_Pos=0;//New_Char;
	unsigned pos;
	int c;
	Start=Start-2;
	
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}

	//Temp_BCB[0]=Branch_Characters[0];Temp_BCB[1]=Branch_Characters[1];Temp_BCB[2]=Branch_Characters[2];Temp_BCB[3]=Branch_Characters[3];
	//memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Tag.Start=M.Forward_Start_LookupX[c];Tag.End=M.Forward_End_LookupX[c];
		Tag.Level=M.Lookupsize + 1; Tag.Skip=0;
		Search_Exact(Current_Tag,Tag,Start,StringLength,revfmi);
	}
	//Branch_Characters[0]=Temp_BCB[0];Branch_Characters[1]=Temp_BCB[1];Branch_Characters[2]=Temp_BCB[2];Branch_Characters[3]=Temp_BCB[3];
	//memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}

void Reverse(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,LEN & L, BWT *fwfmi,MEMX & M)
{	
	unsigned Temp=0;
	int c;
	char New_Char;
	char Mismatch_Count=Tag.Mismatches;
	unsigned pos;
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}
	//Temp_BCR[0]=Branch_Characters[0];Temp_BCR[1]=Branch_Characters[1];Temp_BCR[2]=Branch_Characters[2];Temp_BCR[3]=Branch_Characters[3];
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4) | (Current_Tag[L.STRINGLENGTH-1-3]<<6) | Current_Tag[L.STRINGLENGTH-1-4]<<8 | (Current_Tag[L.STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Tag.Start=M.Backward_Start_LookupX[c];Tag.End=M.Backward_End_LookupX[c];
		Tag.Level=M.Lookupsize + 1;Tag.Skip=0;
		Search_Backwards_Exact(Current_Tag,Tag,L.STRINGLENGTH,L.RH,fwfmi,M);//Backward scan for ?|0
	}
	//Branch_Characters[0]=Temp_BCR[0];Branch_Characters[1]=Temp_BCR[1];Branch_Characters[2]=Temp_BCR[2];Branch_Characters[3]=Temp_BCR[3];
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}
//}------------------------------------------- TWO MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------
unsigned Three_Mismatch(bool batmeth1,char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M,revfmi,fwfmi);//mismatches of the form 0|3, with last mismatch at the end...
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Forward_Pointer=0;

	M.FMIndex=FORWARD;
	if(M.Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End[i],M,revfmi,fwfmi);//Mismatches of the form 3|0, with one mismatch at the first position
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Pointer=0;

	M.FMIndex=REVERSE;
	M.Possible_04_Pointer=M.Mismatches_Forward_Pointer;
	if(M.Mismatches_Forward_Pointer!=M.Possible_03_Pointer)
	{
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=M.Possible_04_Pointer-1;i>=M.Possible_03_Pointer;i--)
			{
				Search_Forwards(M.Mismatches_Forward[i],3,L.LH+1,L.RH,revfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
				if(MAXHITS<=M.Hits) break;
			}
			if(MAXHITS<=M.Hits) return M.Hits;
		}

		Do_Branch=Do_All;*/
		for(int i=M.Possible_04_Pointer-1;i>=M.Possible_03_Pointer;i--)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],3,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=FORWARD;
	M.Possible_40_Pointer=M.Mismatches_Backward_Pointer;
	if(M.Mismatches_Backward_Pointer!=M.Possible_30_Pointer)
	{
		/*if(USEQUALITY)
		{
			Do_Branch=Low_Quality;
			for(int i=Possible_40_Pointer-1;i>=Possible_30_Pointer;i--)
			{
				Search_Backwards(Mismatches_Backward[i],3,LH,LH,fwfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
				if(MAXHITS==Hits) break;
			}
			if(MAXHITS==Hits) return;
		}

		Do_Branch=Do_All;*/
		for(int i=M.Possible_40_Pointer-1;i>=M.Possible_30_Pointer;i--)
		{
			Search_Backwards(batmeth1,Current_Tag,M.Mismatches_Backward[i],3,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}

	//----------------------------------------------------------------------------------------------------------------------------------------
	int c;
	SARange Range;
	M.FMIndex=REVERSE;
	if(M.Lookupsize==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Range.Start=M.Forward_Start_LookupX[c];Range.End=M.Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=M.Lookupsize+1;
	Range.Mismatch_Char=0; Range.Skip=0;
	Search_Forwards_0X(Current_Tag,Range,1,L.LHQL,revfmi);
	Range.Level=1;
	/*if(USEQUALITY)
	{
		Do_Branch=Low_Quality;
		Search_01X(Range,1,LHQL +1,LHQR,revfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
		if(MAXHITS==Hits) return;
	}

	Do_Branch=Do_All;*/
	Search_01X(batmeth1,Current_Tag,Range,1,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M,fwfmi);
	if(MAXHITS<=M.Hits) return M.Hits;
	//----------------------------------------------------------------------------------------------------------------------------------------
	if(M.Lookupsize==3)
	{
		c=Current_Tag[L.LH-1-0] | (Current_Tag[L.LH-1-1]<<2) | (Current_Tag[L.LH-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[L.LH-1-0] | (Current_Tag[L.LH-1-1]<<2) | (Current_Tag[L.LH-1-2]<<4) | (Current_Tag[L.LH-1-3]<<6) | Current_Tag[L.LH-1-4]<<8 | (Current_Tag[L.LH-1-5]<<10);//Use lookup table...
	}

	Range.Start=M.Backward_Start_LookupX[c];Range.End=M.Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=M.Lookupsize+1;//Range.Tag=Actual_Tag;
	Range.Mismatch_Char=0;Range.Skip=0;
	Search_Backwards_Exact_X0( Current_Tag,Range,L.LH,L.LHQR,fwfmi);// ?|0|?
	Range.Level=1;
	/*if(USEQUALITY)
	{
		TRange=Range;
		Do_Branch=Low_Quality;
		Search_10X(TRange,1,LHQL, LHQL,fwfmi);//search for three mismatches of the form 1|2 and stores the candidates for 1|3 
		if(MAXHITS==Hits) return;
	}
	Do_Branch=Do_All;*/
	Search_10X(batmeth1,Current_Tag,Range,1,L.LHQL, L.LHQL,L,MAXHITS,fwfmi,revfmi,M);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	return M.Hits;
}

void Search_10X_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					//Tag.Skip=0;
					Backwards(Current_Tag,Tag,1,L.LH,revfmi,M);
					if(Tag.Start)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=1;
						Search_Forwards(batmeth1,Current_Tag,Tag,3,L.LH+1,L.RH,MAXHITS,revfmi,M,fmi);
					}
					return;

				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else //store mismatches for later use...
		{
			if(M.Right_Mishits_Pointer < M.ARRAY_BOUND)
			{
				Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				M.Right_Mishits[M.Right_Mishits_Pointer]=Tag;
				M.Right_Mishits_Pointer++;
			}
			return;
		}
	}
}
void Search_10X(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L, int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M)
{

	if (!Tag.Start) return;
	struct SARange Range,Temp_Range=Tag;
	int BMHStack_Top=0;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	M.BMHStack[0]=Tag;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_10X_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,fmi,revfmi,M);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								//Temp_Range.Skip=0;
								Backwards(Current_Tag,Temp_Range,1,L.LH,revfmi,M);
								if (Temp_Range.Start)
								{
									Temp_Range.Level=1;
									//Temp_BC[0]=Branch_Characters[0];Temp_BC[1]=Branch_Characters[1];Temp_BC[2]=Branch_Characters[2];Temp_BC[3]=Branch_Characters[3];
									//memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
									Search_Forwards(batmeth1,Current_Tag,Temp_Range,3,L.LH+1,L.RH,MAXHITS,revfmi,M,fmi);
									if(MAXHITS<=M.Hits) return;
									//Branch_Characters[0]=Temp_BC[0];Branch_Characters[1]=Temp_BC[1];Branch_Characters[2]=Temp_BC[2];Branch_Characters[3]=Temp_BC[3];
									//memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
								}

							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						if(M.Right_Mishits_Pointer < M.ARRAY_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							M.Right_Mishits[M.Right_Mishits_Pointer]=Temp_Range;
							M.Right_Mishits_Pointer++;
						}
					}

				} 
			}
		}
	}
	return;
}
void Search_01X(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M,BWT *fwfmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	int FSHStack_Top=0;
	M.FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_01X_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,revfmi,M,fwfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Temp_Range.Level=1;
								//Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
								//memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
								Search_Forwards(batmeth1,Current_Tag,Temp_Range,3,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);
								//memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
								if(MAXHITS<=M.Hits) return;
								//Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						if(M.Left_Mishits_Pointer < M.ARRAY_BOUND)
						{
							Temp_Range.Level += Start+2;
							M.Left_Mishits[M.Left_Mishits_Pointer]=Temp_Range;
							M.Left_Mishits_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}


void Search_01X_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M,BWT *fwfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Forwards(batmeth1,Current_Tag,Tag,3,L.LH+1,L.RH,MAXHITS,fmi,M,fwfmi);
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if(M.Left_Mishits_Pointer < M.ARRAY_BOUND)
			{
				if(!Tag.Skip) Tag.End=Tag.Start;
				Tag.Level += Start+2;
				M.Left_Mishits[M.Left_Mishits_Pointer]=Tag;
				M.Left_Mishits_Pointer++;
			}
			return;
		}
	}
}
//}------------------------------------------- THREE MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------
unsigned Four_Mismatch(bool batmeth1,char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	M.FMIndex=REVERSE;
	SARange Range,TRange;

	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M,revfmi,fwfmi);//mismatches of the form 0|4, with last mismatch at the end...
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.Two_Mismatches_At_End_Forward_Pointer=0;
	M.FMIndex=FORWARD;
	if(M.Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End[i],M,revfmi,fwfmi);//mismatches of the form 0|4, with one mismatch at the start...
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Pointer=0;

	M.FMIndex=REVERSE;
	M.Possible_05_Pointer=M.Mismatches_Forward_Pointer;
	int Wrap=FALSE;
	if(M.Mismatches_Forward_Pointer)
	{
		if (M.Possible_04_Pointer > 80000 ) {M.Mismatches_Forward_Pointer=0;Wrap=TRUE;} //46000
		/*if(USEQUALITY)
		  {
		  Do_Branch=Low_Quality;
		  for(int i=Possible_04_Pointer;i<Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		  {
		  Search_Forwards(Mismatches_Forward[i],4,LH+1,RH,revfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
		  if(MAXHITS==Hits) break;
		  }
		  if(MAXHITS==Hits) continue;
		  }

		  Do_Branch=Do_All;*/
		for(int i=M.Possible_04_Pointer;i<M.Possible_05_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],4,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	if(Wrap) M.Possible_05_Pointer=0;
	M.Mismatches_Forward_Pointer_Last4=M.Mismatches_Forward_Pointer;

	M.FMIndex=FORWARD;
	M.Possible_50_Pointer=M.Mismatches_Backward_Pointer;
	if(M.Mismatches_Backward_Pointer)
	{
		for(int i=M.Possible_50_Pointer-1;i>=M.Possible_40_Pointer;i--)//Mismatches_Backward_Pointer-1;i>=0;i--)
		{
			Search_Backwards(batmeth1,Current_Tag,M.Mismatches_Backward[i],4,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}

	M.FMIndex=REVERSE;
	M.Left_Mishits_Pointer_1=M.Left_Mishits_Pointer;
	if(M.Left_Mishits_Pointer)
	{
		for(int i=0;i<M.Left_Mishits_Pointer;i++)
		{
			if (M.Left_Mishits[i].Level != L.LH+1)
			{
				Search_Exact(Current_Tag,M.Left_Mishits[i],-1,L.LH,revfmi);
				M.Left_Mishits[i].Level=L.LH+1;//search exact does not update level..
			}
			if (M.Left_Mishits[i].Start)
			{
				Search_Forwards(batmeth1,Current_Tag,M.Left_Mishits[i],4,1,L.STRINGLENGTH,MAXHITS,revfmi,M,fwfmi);//find mismatches of the form 022 form, stores possibles of the form 023
			}
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}

	M.Mismatches_Forward_Pointer_Last5=M.Mismatches_Forward_Pointer;
	if( M.Right_Mishits_Pointer)
	{
		for(int i=0;i<M.Right_Mishits_Pointer;i++)
		{

			if(M.Right_Mishits[i].Level!=L.LHQL) 
			{
				Search_Backwards_Exact( Current_Tag,M.Right_Mishits[i],L.LHQL,L.LHQL,fwfmi,M);//finds mismatches of the form 202, stores possibles of the form 203
			}
			if(M.Right_Mishits[i].Start)
			{	
				Backwards(Current_Tag,M.Right_Mishits[i],1,L.LH,revfmi,M);
				if(M.Right_Mishits[i].Start)
				{
					M.Right_Mishits[i].Level=1;
					Search_Forwards(batmeth1,Current_Tag,M.Right_Mishits[i],4,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);
					if(MAXHITS<=M.Hits) break;
				}
			}
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=REVERSE;
	unsigned SOURCELENGTH = fwfmi->textLength;
	int LHQLrx=L.LHQL/2;int c;
	if (L.LHQL % 2) LHQLrx++; int LHQRrx=L.LHQL-LHQLrx;

	if (M.Lookupsize >= LHQLrx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
		}
		Range.Start=M.Forward_Start_LookupX[c];Range.End=M.Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=M.Lookupsize+1;
		Range.Mismatch_Char=0; Range.Skip=0;
	}
	Search_Forwards_0X(Current_Tag,Range,1,LHQLrx,revfmi);
	Range.Level=1;
	Search_01LX(batmeth1,Current_Tag,Range,1,LHQLrx +1,LHQRrx,MAXHITS,revfmi,L,M,fwfmi);
	if(MAXHITS<=M.Hits) return M.Hits;

//--------------------------------
	if (M.Lookupsize >= LHQRrx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[L.LHQL-1-0] | (Current_Tag[L.LHQL-1-1]<<2) | (Current_Tag[L.LHQL-1-2]<<4);// | (Current_Tag[LH-1-3]<<6) | Current_Tag[LH-1-4]<<8 | (Current_Tag[LH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[L.LHQL-1-0] | (Current_Tag[L.LHQL-1-1]<<2) | (Current_Tag[L.LHQL-1-2]<<4) | (Current_Tag[L.LHQL-1-3]<<6) | Current_Tag[L.LHQL-1-4]<<8 | (Current_Tag[L.LHQL-1-5]<<10);//Use lookup table...
		}

		Range.Start=M.Backward_Start_LookupX[c];Range.End=M.Backward_End_LookupX[c];Range.Mismatches=0;Range.Level=M.Lookupsize+1;//Range.Tag=Actual_Tag;
		Range.Mismatch_Char=0;Range.Skip=0;
	}
	Search_Backwards_Exact_X0( Current_Tag,Range,L.LHQL,LHQRrx,fwfmi);// ?|0|?
	Range.Level=1;
	//Do_Branch=Do_All;
	Search_10LX(batmeth1,Current_Tag,Range,1,LHQLrx, LHQLrx,MAXHITS,fwfmi,revfmi,L,M);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	return M.Hits;
}

void Search_10LX_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fwfmi->inverseSa0) Index--;//adjust for missing $
		Now=fwfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		//if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fwfmi->cumulativeFreq[Now] + BWTOccValue(fwfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					//Tag.Skip=0;
					Backwards(Current_Tag,Tag,1,L.LHQL,revfmi,M);
					if(Tag.Start)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=1;
						Search_Half_Tag_11X(batmeth1,Current_Tag,Tag,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M,fwfmi);
					}
					return;

				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else 
		{
			return;
		}
	}
}

void Search_10LX(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
{
	
	if (!Tag.Start) return;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	struct SARange Range,Temp_Range=Tag;
	int BMHStack_Top=0;
	M.BMHStack[0]=Tag;
	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_10LX_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,fwfmi,revfmi,L,M);
			//if (M.Last_Mismatch_Written ==8 ) return;
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fwfmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start , Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								//Temp_Range.Skip=0;
								Backwards(Current_Tag,Temp_Range,1,L.LHQL,revfmi,M);
								if (Temp_Range.Start)
								{
									Temp_Range.Level=1;
									Search_Half_Tag_11X(batmeth1,Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M,fwfmi);
									//if (M.Last_Mismatch_Written ==8 ) return;
									if(MAXHITS<=M.Hits) return;
								}

							}
							else continue;
						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	return;
}
void Search_01LX_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,LEN & L,MEMX & M,BWT *fwfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=1;
					Search_Half_Tag_11X(batmeth1,Current_Tag,Tag,2,L.LHQL +1,L.LHQR,L,MAXHITS,fmi,M,fwfmi);
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			return;
		}
	}
}

void Search_01LX(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *revfmi,LEN & L,MEMX & M,BWT *fwfmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	M.FSHStackX0X[0]=Tag;
	struct SARange Range,Temp_Range;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];

	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStackX0X[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_01LX_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,revfmi,L,M,fwfmi);
			//if (M.Last_Mismatch_Written ==8 ) return;
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								Temp_Range.Level=1;
								Search_Half_Tag_11X(batmeth1,Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M,fwfmi);
								//if (M.Last_Mismatch_Written ==8 ) return;
								if(MAXHITS<=M.Hits) return;
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSHStackX0X[FSHStack_Top]=Temp_Range;
						}
					}
				} 
			}

		}
	}
	return;
}


void Search_11X(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M,BWT *fwfmi)
{

	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	M.FSHStackX0X[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStackX0X[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);
		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;
				Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
				Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				if (Current_Tag[Temp_Range.Level+Start]!=Branch)
				{
					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
					Temp_Range.Mismatches++;
				}

				if (Temp_Range.Mismatches<=1)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength)
					{
						Temp_Range.Level=1;
						//Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
						//memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));
						Search_Half_Tag_11X(batmeth1,Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M,fwfmi);
						if(MAXHITS<=M.Hits) return;
						//memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
						//Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];
					}
					else
					{
						FSHStack_Top++;//Push range
						Temp_Range.Level++;
						M.FSHStackX0X[FSHStack_Top]=Temp_Range;
					}
				}
			} 
		}

	}
	return;
}

void Search_Half_Tag_11X_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M,BWT *fwfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)// && Tag.Mismatches==Count)
			{
				if (Tag.Mismatches != Count) return;
				if(!Tag.Skip) Tag.End=Tag.Start;
				Tag.Level=1;
				Search_Forwards(batmeth1,Current_Tag,Tag,4,L.LH+1,L.RH,MAXHITS,fmi,M,fwfmi);
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if (M.Left_Mishits_Pointer < M.ARRAY_BOUND)//Bounds Check...
			{
				if(!Tag.Skip) Tag.End=Tag.Start;
				//Tag.Level=LHQR+Tag.Level;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				//Tag.Level=LHQR+Tag.Level;//+1;
				Tag.Level += Start+2;
				M.Left_Mishits[M.Left_Mishits_Pointer]=Tag;
				M.Left_Mishits_Pointer++;
			}
			return;
		}
	}
}

void Search_Half_Tag_11X(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M,BWT *fwfmi)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	int FSHStack_Top=0;
	M.FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Half_Tag_11X_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,revfmi,M,fwfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength )//&& Temp_Range.Mismatches==Count)
						{
							if (Temp_Range.Mismatches != Count) continue;
							Temp_Range.Level=1;
							//Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
							//memcpy(Temp_Branch_Ranges,Branch_Ranges,4*sizeof(SARange));
							Search_Forwards(batmeth1,Current_Tag,Temp_Range,4,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);
							if(MAXHITS<=M.Hits) return;
							//memcpy(Branch_Ranges,Temp_Branch_Ranges,4*sizeof(SARange));
							//Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else //store mismatches for later use...
					{
						if(M.Left_Mishits_Pointer < M.ARRAY_BOUND) 
						{
							//Temp_Range.Level=LHQR+Temp_Range.Level;
							//if (Temp_Range.Level!=StringLength) Temp_Range.Level++;
							//Temp_Range.Level=LHQR+Temp_Range.Level;//+1;
							Temp_Range.Level += Start+2;
							M.Left_Mishits[M.Left_Mishits_Pointer]=Temp_Range;
							M.Left_Mishits_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}
//}------------------------------------------- FOUR MISMATCH ---------------------------------------------------------------------------------------------

//{------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------
unsigned Five_Mismatch(bool batmeth1,char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{

	SARange Range,TRange;
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M,revfmi,fwfmi);//mismatches of the form 0|4, with last mismatch at the end...
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Forward_Pointer=0;


	M.FMIndex=FORWARD;
	if(M.Two_Mismatches_At_End_Pointer)
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Pointer;i++)
		{
			if(batmeth1) Print_LocationX(M.Two_Mismatches_At_End[i],M);
			else Print_LocationX(M.Two_Mismatches_At_End[i],M,revfmi,fwfmi);//mismatches of the form 0|5, with one mismatch at the start...
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	M.Two_Mismatches_At_End_Pointer=0;

	M.FMIndex=REVERSE;
	if(M.Mismatches_Forward_Pointer)
	{
		for(int i=M.Possible_05_Pointer;i<M.Mismatches_Forward_Pointer_Last4;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],5,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for possible four mismatches of the form 0|5
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}


	M.FMIndex=REVERSE;
	if(M.Mismatches_Forward_Pointer)
	{
		for(int i=M.Mismatches_Forward_Pointer_Last4;i<M.Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],5,1,L.STRINGLENGTH,MAXHITS,revfmi,M,fwfmi);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=REVERSE;
	if(M.Mismatches_Forward_Pointer)
	{
		for(int i=M.Mismatches_Forward_Pointer_Last5;i<M.Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(batmeth1,Current_Tag,M.Mismatches_Forward[i],5,L.LH+1,L.RH,MAXHITS,revfmi,M,fwfmi);//scan for possible four mismatches of the form 0|5
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=FORWARD;
	if(M.Mismatches_Backward_Pointer!=M.Possible_50_Pointer)
	{
		for(int i=M.Mismatches_Backward_Pointer-1;i>=M.Possible_50_Pointer;i--)
		{
			Search_Backwards(batmeth1,Current_Tag,M.Mismatches_Backward[i],5,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}
	M.FMIndex=FORWARD;
	if (M.Possible_20_Pointer)
	{
		for(int i=M.Possible_20_Pointer-1;i>=0;i--)
		{
			if(M.Possible_20[i].Level!=L.RHQL) 
			{
				M.Possible_20[i].Level++; 
				Search_Backwards_Exact( Current_Tag,M.Possible_20[i],L.LH+L.RHQL,L.RHQL,fwfmi,M);
			}

			if(M.Possible_20[i].Start)
			{	
				M.Possible_20[i].Level=1;
				Search_Backwards(batmeth1,Current_Tag,M.Possible_20[i],5,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);
				if(MAXHITS<=M.Hits) break;
			}
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}

	if(M.Possible_02_Pointer)
	{
		for(int i=0;i<M.Possible_02_Pointer;i++)
		{

			if(M.Possible_02[i].Level!=L.RHQR) 
			{
				Search_Exact(Current_Tag,M.Possible_02[i],L.LH + L.RHQL-1 ,L.RHQR,revfmi);//finds mismatches of the form 202, stores possibles of the form 203
			}
			if(M.Possible_02[i].Start) 
			{
				Reverse(Current_Tag,M.Possible_02[i],L.STRINGLENGTH,L.RH,L,fwfmi,M);
				if(M.Possible_02[i].Start)
				{
					M.Possible_02[i].Level=1;
					Search_Backwards(batmeth1,Current_Tag,M.Possible_02[i],5,L.LH,L.LH,MAXHITS,fwfmi,M,revfmi);
				}
			}
			if(MAXHITS<=M.Hits) break;
		}

		if(MAXHITS<=M.Hits) return M.Hits;
	}

	int RHQLlx=L.RHQR/2;
	unsigned SOURCELENGTH = fwfmi->textLength;int c;
	if (L.RHQR % 2) RHQLlx++; int RHQRlx=L.RHQR-RHQLlx;
	if (M.Lookupsize >= RHQRlx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[L.STRINGLENGTH-1-0] | (Current_Tag[L.STRINGLENGTH-1-1]<<2) | (Current_Tag[L.STRINGLENGTH-1-2]<<4) | (Current_Tag[L.STRINGLENGTH-1-3]<<6) | Current_Tag[L.STRINGLENGTH-1-4]<<8 | (Current_Tag[L.STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		Range.Start=M.Backward_Start_LookupX[c];Range.End=M.Backward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=M.Lookupsize+1;//Range.Tag=Actual_Tag;
		Range.Mismatch_Char=0;Range.Skip=0;
	}

	Search_Backwards_Exact_X0( Current_Tag,Range,L.STRINGLENGTH,RHQRlx,fwfmi);// ?|?|0
	Range.Level=1;

	Search_Backwards_XL10(batmeth1,Current_Tag,Range,1,L.LH + L.RHQL+RHQLlx, RHQLlx,MAXHITS,fwfmi,L,M,revfmi);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	if(MAXHITS<=M.Hits) return M.Hits;


	if (M.Lookupsize >= RHQLlx)
	{
		Range.Start=1;Range.End=SOURCELENGTH;Range.Mismatches=0;Range.Level=1;
		Range.Mismatch_Char=0;Range.Skip=0;
	}
	else
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[L.LH+L.RHQL+0] | (Current_Tag[L.LH+L.RHQL+1]<<2) | (Current_Tag[L.LH+L.RHQL+2]<<4);// | (Current_Tag[LH+3]<<6) | Current_Tag[LH+4]<<8 | (Current_Tag[LH+5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[L.LH+L.RHQL+0] | (Current_Tag[L.LH+L.RHQL+1]<<2) | (Current_Tag[L.LH+L.RHQL+2]<<4) | (Current_Tag[L.LH+L.RHQL+3]<<6) | Current_Tag[L.LH+L.RHQL+4]<<8 | (Current_Tag[L.LH+L.RHQL+5]<<10);//Use lookup table...
		}
		Range.Start=M.Forward_Start_LookupX[c];Range.End=M.Forward_End_LookupX[c];
		Range.Mismatches=0;Range.Level=M.Lookupsize+1; Range.Skip=0;
		Range.Mismatch_Char=0;
	}

	Search_Forwards_0X(Current_Tag,Range,L.LH+L.RHQL+1,RHQLlx,revfmi);
	Range.Level=1;TRange=Range;
	Search_XL01(batmeth1,Current_Tag,TRange,1,L.LH + L.RHQL+RHQLlx +1,RHQRlx,MAXHITS,fwfmi,revfmi,L,M);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	return M.Hits;
}

void Search_XL01_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* fwfmi,BWT *revfmi,LEN & L,MEMX & M)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= revfmi->inverseSa0) Index--;//adjust for missing $
		Now=revfmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		//if ( !Do_Branch[Start+Tag.Level] && Now != Current_Tag[Start+Tag.Level]) return; //do not bend these nuces...
		Tag.Start = revfmi->cumulativeFreq[Now] + BWTOccValue(revfmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start+Tag.Level] != Now) 
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
		}
		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)
				{
					int T=L.RH;L.RH=L.RHQR;
					Reverse(Current_Tag,Tag,L.STRINGLENGTH,L.RH,L,fwfmi,M);
					L.RH=T;
					if(Tag.Start)
					{
						if(!Tag.Skip) Tag.End=Tag.Start;
						Tag.Level=L.RHQR+1;//Temp_Range.Level=1;
						Search_Half_Tag_X11(batmeth1,Current_Tag,Tag,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fwfmi,M,revfmi);
					}
					return;
				}
				else return;
			}
			else {Tag.Level++;continue;}
		} 
		else
		{
			if(M.Possible_02_Pointer < M.END_BOUND)
			{
				if(!Tag.Skip)Tag.End=Tag.Start;
				if (Tag.Level!=StringLength) Tag.Level++; 
				M.Possible_02[M.Possible_02_Pointer]=Tag;
				M.Possible_02_Pointer++;
			}
			return;
		}
	}
}
void Search_XL01(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSHStack_Top=0;
	M.FSHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];

	while(FSHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSHStack[FSHStack_Top];
		FSHStack_Top--;		//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_XL01_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,fwfmi,revfmi,L,M);
			//if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			//if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)
							{
								int T=L.RH;L.RH=L.RHQR;
								Reverse(Current_Tag,Temp_Range,L.STRINGLENGTH,L.RH,L,fwfmi,M);
								L.RH=T;
								if(Temp_Range.Start)
								{
									Temp_Range.Level=L.RHQR+1;//Temp_Range.Level=1;
									Search_Half_Tag_X11(batmeth1,Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fwfmi,M,revfmi);
									//if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
									if(MAXHITS<=M.Hits) return;

								}
							}
							else continue;
						}
						else
						{
							FSHStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSHStack[FSHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(M.Possible_02_Pointer < M.END_BOUND)
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							M.Possible_02[M.Possible_02_Pointer]=Temp_Range;
							M.Possible_02_Pointer++;
						}
					}
				} 
			}

		}
	}
	return;
}
void Search_Backwards_XL10_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,LEN & L, MEMX & M,BWT *revfmi)
{

	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if ( !Do_Branch[Start-Tag.Level] && Now != Current_Tag[Start-Tag.Level]) return; //do not bend these nuces...
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)//only one mismatch allowed here...
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches)//a tag of the form ?|1|0 , remove zero mismatch
				{
					if(!Tag.Skip) Tag.End=Tag.Start;
					Tag.Level=L.RHQR+1;//Temp_Range.Level=1;
					Search_Half_Tag_X11(batmeth1,Current_Tag,Tag,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M,revfmi);
					return;
				}
				else return;
			}
			else { Tag.Level++;continue; }
		} 
		else
		{
			if(M.Possible_20_Pointer < M.END_BOUND)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				//if (Tag.Level!=StringLength) Tag.Level++; 
				M.Possible_20[M.Possible_20_Pointer]=Tag;
				M.Possible_20_Pointer++;
			}
			return;
		}
	}
}
void Search_Backwards_XL10(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *fmi,LEN & L,MEMX & M,BWT *revfmi)
{
	if (!Tag.Start) return;
	int BMHStack_Top=0;
	M.BMHStack[0]=Tag;
	struct SARange Range,Temp_Range;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];

	while(BMHStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMHStack[BMHStack_Top];
		BMHStack_Top--;	//Pop the range

		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_XL10_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,MAXHITS,fmi,L,M,revfmi);
			//if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			//if (SEED && Range.End-Range.Start >20) continue;//@@@@
			Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;
					Temp_Range.End = Branch_Ranges[Branch].End;

					if (Current_Tag[Start-Temp_Range.Level] != Branch)//only one mismatch allowed here...
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches)//a tag of the form ?|1|0
							{
								Temp_Range.Level=L.RHQR+1;//Temp_Range.Level=1;
								Search_Half_Tag_X11(batmeth1,Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M,revfmi);
								//if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
								if(MAXHITS<=M.Hits) return;
							}
							else continue;

						}
						else
						{
							BMHStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMHStack[BMHStack_Top]=Temp_Range;
						}
					}
					else
					{
						if (M.Possible_20_Pointer < M.END_BOUND)
						{
							M.Possible_20[M.Possible_20_Pointer]=Temp_Range;
							M.Possible_20_Pointer++;
						}
					}
				} 
			}
		}
	}
	return;
}

void Search_X11(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,unsigned MAXHITS, LEN & L,BWT *fmi,MEMX & M,BWT *revfmi)
{
	if (!Tag.Start) return;
	int BMStack_Top=0;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	M.BMStack_X11[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMStack_X11[BMStack_Top];
		BMStack_Top--;	//Pop the range
		Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);

		for(int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])//This character actually branches
			{
				Temp_Range=Range;//adjust
				Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
				Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

				if (Current_Tag[Start-Temp_Range.Level] != Branch)
				{
					Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
					Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
					Temp_Range.Mismatches++;
				}

				if (Temp_Range.Mismatches<=1)//we are guaranteed a valid SA range, check only for mismatches
				{
					if(Temp_Range.Level== StringLength)
					{
						if(Temp_Range.Mismatches)
						{
							Temp_Range.Level++;
							//Temp_BC2[0]=Branch_Characters[0];Temp_BC2[1]=Branch_Characters[1];Temp_BC2[2]=Branch_Characters[2];Temp_BC2[3]=Branch_Characters[3];
							//memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));//X|8
							Search_Half_Tag_X11(batmeth1,Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M,revfmi);
							if(MAXHITS<=M.Hits) return;
							//memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							//Branch_Characters[0]=Temp_BC2[0];Branch_Characters[1]=Temp_BC2[1];Branch_Characters[2]=Temp_BC2[2];Branch_Characters[3]=Temp_BC2[3];

						}
						else continue;
					}
					else
					{
						BMStack_Top++;//Push range
						Temp_Range.Level++;
						M.BMStack_X11[BMStack_Top]=Temp_Range;
					}
				}
			} 
		}
	}
	return;
}

void Search_Half_Tag_X11_OneSA(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M,BWT *revfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}


	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		//if (!Do_Branch[Start-Tag.Level] && Current_Tag[Start-Tag.Level]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength && Tag.Mismatches==2)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				Tag.Level=1;
				Search_Backwards(batmeth1,Current_Tag,Tag,5,L.LH,L.LH,MAXHITS,fmi,M,revfmi);//LH,fwfmi);
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else return;
	}
}

void Search_Half_Tag_X11(bool batmeth1,char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M,BWT *revfmi)
{
	if (!Tag.Start) return;

	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	int BMStack_Top=0;
	M.BMStack_X11H[0]=Tag;
	struct SARange Range,Temp_Range;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMStack_X11H[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Half_Tag_X11_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,fmi,M,revfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start-Temp_Range.Level);
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=2)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength && Temp_Range.Mismatches == Count)
						{
							Temp_Range.Level=1;
							//Temp_BC1[0]=Branch_Characters[0];Temp_BC1[1]=Branch_Characters[1];Temp_BC1[2]=Branch_Characters[2];Temp_BC1[3]=Branch_Characters[3];
							//memcpy(Temp_Branch_Ranges2,Branch_Ranges,4*sizeof(SARange));//X|16
							Search_Backwards(batmeth1,Current_Tag,Temp_Range,5,L.LH,L.LH,MAXHITS,fmi,M,revfmi);//LH,fwfmi);
							if(MAXHITS<=M.Hits) return;
							//memcpy(Branch_Ranges,Temp_Branch_Ranges2,4*sizeof(SARange));
							//Branch_Characters[0]=Temp_BC1[0];Branch_Characters[1]=Temp_BC1[1];Branch_Characters[2]=Temp_BC1[2];Branch_Characters[3]=Temp_BC1[3];
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMStack_X11H[BMStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	return;
}
//}------------------------------------------- FIVE MISMATCH ---------------------------------------------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Guess_Orientation
 *  Description:  Guesses the best direction. If hits found, returns the hit count. Store
 *  		  guessed direction in G.
 * =====================================================================================
 */
unsigned Guess_Orientation(bool batmeth1,BWT* fwfmi,BWT* revfmi,MEMX & MF,MEMX & MC,LEN & L,GUESS & G,BATREAD & B)
{
	int Complement_Length;
	int Org_Length;
	SARange Range;
	MF.FMIndex=REVERSE;MC.FMIndex=REVERSE;
	unsigned Hits=0;
	int c;

//-------------------------------------------------- NORMAL SCAN --------------------------------------------------------------------
	//Exact_Match_Forward=MF.Exact_Match_Forward;
	//Current_Tag=Original+IGNOREHEAD;
	char *Current_Tag=B.Forward;
	MF.Guessed=FALSE;MC.Guessed=FALSE;//was the particular strand scanned during guessing?
	MF.Strand='+';
	if(MF.Lookupsize ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}

	Range.Start=MF.Forward_Start_LookupX[c];Range.End=MF.Forward_End_LookupX[c];
	Range.Mismatches=0;Range.Level=MF.Lookupsize+1;Range.Skip=0;Range.Mismatch_Char=0;
	Search_Forwards_Exact(batmeth1,Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,MF,L.LH,revfmi,fwfmi);//Find exact matches and report... if not found get the range for 0|?
	MF.Guessed=TRUE;
	if(MF.Hits) 
	{
		Hits += MF.Hits;
		G.Guessed=MF;//Original+IGNOREHEAD;
		G.Guessed_Read=B.Forward;
		//G.Guessed_Quality=Low_QualityF+IGNOREHEAD;
		G.Guess_Complement=MC;
		G.Guessed_ReadC=B.Complement;
		//G.Guessed_Complement_Quality=Low_QualityC+IGNOREHEAD;
		//if (MAX_MISMATCHES != 0) {Hit_In=1;return FALSE;}
		//if (!ROLLOVER ) return FALSE;
		return Hits;
	}
	Org_Length=Range.Level;
//--------------------------------------------------REVERSE COMPLEMENT SCAN --------------------------------------------------------------------
	//Exact_Match_Forward=MC.Exact_Match_Forward;
	Current_Tag=B.Complement;
	MC.Strand='-';
	//NLocations=NLocationsC;

	if(MC.Lookupsize ==3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}

	Range.Start=MC.Forward_Start_LookupX[c];Range.End=MC.Forward_End_LookupX[c];Range.Mismatches=0;Range.Level=MC.Lookupsize+1;
	Range.Mismatch_Char=0;Range.Skip=0;
	Search_Forwards_Exact(batmeth1,Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,MC,L.LH,revfmi,fwfmi);//Find exact matches and report... if not found get the range for 0|?
	Hits += MC.Hits;
	Complement_Length=Range.Level;
	MC.Guessed=TRUE;
//-------------------------------------------------- EXACT SCAN END --------------------------------------------------------------------
	if (Complement_Length > Org_Length)//first try the complement
	{
		//Max_Seeked=Complement_Length;
		//Delta=Complement_Length-Org_Length;
		G.Guessed=MC;//Complement;
		G.Guessed_Read=B.Complement;
		//Guessed_NLocation=NLocationsC;
		//Guessed_Quality=Low_QualityC+IGNOREHEAD;
		G.Guess_Complement=MF;//Original+IGNOREHEAD;
		G.Guessed_ReadC=B.Forward;
		//Guessed_Complement_Quality=Low_QualityF+IGNOREHEAD;
	}
	else
	{

		//Delta=Org_Length-Complement_Length;
		G.Guessed=MF;
		G.Guessed_Read=B.Forward;
		//Guessed_Quality=Low_QualityF+IGNOREHEAD;
		G.Guess_Complement=MC;
		G.Guessed_ReadC=B.Complement;
		//Guessed_NLocation_Complement=NLocationsC;
	}
	return Hits;

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Forwards_Exact
 *  Description:  search for an exact match of the hit, from Current_Tag[Start+Tag.Level ..StringLength]
 *  		  BWT is revfmi
 * =====================================================================================
 */
void Search_Forwards_Exact(bool batmeth1,char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi,MEMX & M,int LH,BWT *revfmi,BWT *fwfmi)
{
	int Level;
	unsigned Branch_Characters[4];
	unsigned long Index,Now,First,Last;

	M.Exact_Match_Forward[Start+LH].Start=0;
	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}

			for(;;)
			{
				
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start+Tag.Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}
					M.Exact_Match_Forward[Start+Tag.Level]=Tag;
					if (!Tag.Skip) M.Exact_Match_Forward[Start+Tag.Level].End=Tag.Start;

					if(Tag.Level== StringLength)
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						if(batmeth1) Print_LocationX(Tag,M);
						else Print_LocationX(Tag,M,revfmi,fwfmi);
						return;	
					}
					else {Tag.Level++;continue;}
				} 
				else//mismatch...
				{
					M.Exact_Match_Forward[0]=Tag;//save old location for heuristics...
					Tag.Start=0;
					M.Exact_Match_Forward[Start+Tag.Level]=Tag;
					Tag.Level--;
					return;	
				}
			}
		}
		else//SA range has sevaral possible hits... 
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Tag.Level+Start];
				if (Branch_Characters[Now])//we have a match... 
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				}
				else//mismatch..
				{
					M.Exact_Match_Forward[0]=Tag;//save old location for heuristics...
					Tag.Start=0;
				}
			} 
			else
			{
				M.Exact_Match_Forward[0]=Tag;//save old location for heuristics...
				Get_SARange_Fast(Current_Tag[Start+Tag.Level],Tag,fmi);
			}

			M.Exact_Match_Forward[Tag.Level+Start]=Tag;
			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					if(batmeth1) Print_LocationX(Tag,M);
					else Print_LocationX(Tag,M,revfmi,fwfmi);
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else//Mismatch
			{
				M.Exact_Match_Forward[Start+Tag.Level]=Tag;
				Tag.Level--;
				return;
			}

		}
	}
}

void Search_Forwards(bool batmeth1,const char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* & revfmi,MEMX & M,BWT *fwfmi)
{
	unsigned Branch_Characters[4];
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	int FSStack_Top=0;
	M.FSSStack[0]=Tag;
	SARange Branch_Ranges[4];
	struct SARange Range,Temp_Range;

	while(FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSSStack[FSStack_Top];
		FSStack_Top--;		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Forwards_OneSA(batmeth1,Current_Tag,Range,Count,Start,StringLength,revfmi,M,revfmi,fwfmi);
			if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches == Count) //dont print exact matches
							{
								if(batmeth1) Print_LocationX(Temp_Range,M);
								else Print_LocationX(Temp_Range,M,revfmi,fwfmi);
								if (MAXHITS<=M.Hits) return;
							}
							else continue;
						}
						else 
						{

							FSStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSSStack[FSStack_Top]=Temp_Range;
						}
					}
					else
					{
						if(5 >Count)//store only for one mismatch... and last node will not branch
						{
							if (Temp_Range.Level!=StringLength) Temp_Range.Level++; 
							else //2 mismatches with the last at the end...
							{
								if(M.Two_Mismatches_At_End_Forward_Pointer < M.END_BOUND)
								{
									M.Two_Mismatches_At_End_Forward[M.Two_Mismatches_At_End_Forward_Pointer]=Temp_Range;
									M.Two_Mismatches_At_End_Forward_Pointer++;
								}
								continue;
							}
							if(M.Mismatches_Forward_Pointer < M.ARRAY_BOUND)
							{
								M.Mismatches_Forward[M.Mismatches_Forward_Pointer]=Temp_Range;
								M.Mismatches_Forward_Pointer++;
							}
						}
						continue;
					}
				} 
			}
		}
	}
	return;
}

void Search_Forwards_OneSA(bool batmeth1,const char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M,BWT *revfmi,BWT *fwfmi)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		//if (!Do_Branch[Tag.Level+Start] && Current_Tag[Tag.Level+Start]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)//avoid printing exact matches...
				{
					if (!Tag.Skip) Tag.End=Tag.Start;
					if(batmeth1) Print_LocationX(Tag,M);
					else Print_LocationX(Tag,M,revfmi,fwfmi);
				}
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else//log 2 mismatches 
		{
			if(5 > Count)//store only for one mismatch... later report these on a seperate stack..
			{
				if (!Tag.Skip) Tag.End=Tag.Start;//possibly two mismatch exists..
				if (Tag.Level != StringLength) Tag.Level++; 
				else //2 mismatches occuring in last position...
				{
					if(M.Two_Mismatches_At_End_Forward_Pointer < M.END_BOUND)
					{
						//if(Tag.Skip) Tag.Start=Tag.End;
						M.Two_Mismatches_At_End_Forward[M.Two_Mismatches_At_End_Forward_Pointer]=Tag;
						M.Two_Mismatches_At_End_Forward_Pointer++;
					}
					return;
				}
				if(M.Mismatches_Forward_Pointer < M.ARRAY_BOUND)
				{
					M.Mismatches_Forward[M.Mismatches_Forward_Pointer]=Tag;
					M.Mismatches_Forward_Pointer++;
				}
			}
			return;
		}
	}
}

void Branch_Detect (const struct SARange Tag,BWT *fmi,int Start,unsigned* Branch_Characters,SARange* Branch_Ranges)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/ if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
			/*if ( !Do_Branch[Tag.Level+Start] && Branch != Current_Tag[Start+Tag.Level]) 
			{
				Branch_Characters[Branch]=0; //do not bend these nuces...
			}
			else*/
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
			}
		}

	}
}


void Search_Backwards_Exact(const char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi,MEMX & M)
{

	int Level;
	unsigned long Index,Now,First,Last;
	if (!Tag.Start) return;
	unsigned Branch_Characters[4];

	for(;;)	
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			Level=Tag.Level;
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
			{
				Tag.Skip++;
				Tag.End=Tag.Start;
			}
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start-Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0) 
					{
						Tag.Skip++;Tag.End=Tag.Start;
					}

					//Exact_Match_Backward[Level]=Tag;
					//if (!Tag.Skip) Exact_Match_Backward[Level].End=Tag.Start;
					if(Level== StringLength)//no need to print as we are still halfway..
					{
						if (!Tag.Skip) Tag.End=Tag.Start;
						return;	
					}
					else {Level++;continue;}
				} 
				else
				{
					Tag.Start=0;//Tag.End=0;
					return;	
				}
			}
		}
		else 
		{
			M.Exact_Match_Backward[Tag.Level]=Tag;
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 
				for (unsigned Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;	
				}

				Now=Current_Tag[Start-Tag.Level];
				if (Branch_Characters[Now])//we have a match... 
				{

					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges
				
				}
				else
				{
					Tag.Start=0;//Tag.End=0;
				}
			}
			else Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fmi);

			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					return;
				}
				else {Tag.Level++;continue;}
			} 
			else
			{
				return;
			}

		}
	}
}
//}-----------------------------------  SCAN ROUTINES ---------------------------------------------------------

//{-----------------------------------  PRINT ROUTINES ---------------------------------------------------------
void Write_Header(HEADER & Header,FILE* Output_File,BATPARAMETERS P,INFILE F,LEN L)
{
#define GISMODE 100
	int SOLIDMARK;
	Header.ID[0]='B';Header.ID[1]='A';Header.ID[2]='T';
	Header.MAXHITS=P.MAXHITS;
	Header.FILETYPE=F.FILETYPE;
	Header.HITMODE = P.Patternfile_Count ? PAIRED_END: SINGLE_END;//HITMODE;
	Header.IGNOREHEAD = L.IGNOREHEAD;
	Header.Tag_Length=L.STRINGLENGTH;
	if (F.SOLID) SOLIDMARK=100; else SOLIDMARK=0;
	Header.Index_Count=P.ONEFMINDEX+SOLIDMARK;
	//if (GIS) PRINT_DESC=GISMODE;
	Header.Print_Desc=GISMODE;
	fwriteX(&Header,sizeof(Header),1,Output_File);
	fwriteX(&P.MAX_MISMATCHES,sizeof(char),1,Output_File);
	fwriteX(&F.TAG_COPY_LEN,sizeof(int),1,Output_File);
	fwriteX(&P.ROLLOVER,sizeof(char),1,Output_File);
	fwriteX(&P.SCANBOTH,sizeof(char),1,Output_File);
}
void Print_LocationX(SARange & Tag, MEMX & M)
{
	unsigned Gap;
	int In_Mismatch= Tag.Mismatches<5 ? Tag.Mismatches : 6;
	Tag.FMIndex=M.FMIndex;
	Tag.Strand=M.Strand;
	M.Hit_Array[M.Hit_Array_Ptr++]=Tag;
	if (M.Hit_Array_Ptr >= 999) M.Hit_Array_Ptr=999;
	if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
	Gap++;M.Hits += Gap;
	M.Stats[In_Mismatch] +=Gap;
}
void Print_LocationX(SARange & Tag, MEMX & M,BWT *revfmi,BWT *fwfmi)
{
	unsigned Gap;
	int In_Mismatch= Tag.Mismatches<5 ? Tag.Mismatches : 6;
	Tag.FMIndex=M.FMIndex;
	Tag.Strand=M.Strand;

	assert (Tag.Skip || Tag.End>=Tag.Start);

	/*if(!M.Extend || M.Extend==PLUS)
	{
		if ( FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse(Tag,M.Extend?(M.L.STRINGLENGTH/2):(M.L.STRINGLENGTH),revfmi,M);
		}
	}
	else */if(M.Extend==PLUS)
	{
		if ( FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse(Tag,M.Extend?(M.L.STRINGLENGTH/2):(M.L.STRINGLENGTH),revfmi,M);
		}
		Tag.Level=M.L.STRINGLENGTH/2+1;
		Tag=Seed_ExtendF(M.Current_Tag,Tag,10,1,M.L.STRINGLENGTH,100,revfmi,M);
		if(Tag.Start)
		{
			if ( FORWARD==Tag.FMIndex )// forward search index...
			{
				Convert_To_Reverse(Tag,(M.L.STRINGLENGTH),revfmi,M);
			}
			assert(Tag.Start);
		}
	}
	else if(M.Extend==MINUS)
	{
		//if ( FORWARD!=Tag.FMIndex )// forward search index...
		{
			Convert_To_Fwd(M.Current_Tag,Tag,M.L.STRINGLENGTH/2,fwfmi,M);
			Tag.FMIndex=FORWARD;
		}
		if(Tag.Mismatches)
		{
			for( int i=0;i<Tag.Mismatches;i++)
			{
				Tag.Mismatch_Pos[i]=Tag.Mismatch_Pos[i]+M.L.STRINGLENGTH/2;
			}

		}
		Tag.Level=M.L.STRINGLENGTH/2+1;
		Tag=Seed_ExtendB(M.Current_Tag+M.L.STRINGLENGTH/2,Tag,10,M.L.STRINGLENGTH,M.L.STRINGLENGTH,100,fwfmi,M);
		if(Tag.Start)
		{
			char Temp_CS[MAXDES];

			for (int i=0;i<M.L.STRINGLENGTH/2;i++) {Temp_CS[i]=M.Current_Tag[i];M.Current_Tag[i]=M.Current_Tag[i+M.L.STRINGLENGTH/2];}
			for (int i=0;i<M.L.STRINGLENGTH/2;i++) {M.Current_Tag[i+M.L.STRINGLENGTH/2]=Temp_CS[i];}
			if(Tag.Mismatches)
			{
				for( int i=0;i<Tag.Mismatches;i++)
				{
					Tag.Mismatch_Pos[i]=Tag.Mismatch_Pos[i]+M.L.STRINGLENGTH/1+1;
				}

			}

			if ( FORWARD==Tag.FMIndex )// forward search index...
			{
				Convert_To_Reverse(Tag,(M.L.STRINGLENGTH),revfmi,M);
			}
			assert(Tag.Start);
		}
	}

	assert (Tag.Skip || Tag.End>=Tag.Start);
	if (Tag.Start)
	{
		unsigned Conversion_Factor=revfmi->textLength-M.L.STRINGLENGTH;
		if (Tag.Skip) Tag.Start=Tag.End;
		if ( FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse(Tag,M.Extend?(M.L.STRINGLENGTH/2):(M.L.STRINGLENGTH),revfmi,M);
		}
		if (M.Hit_Array_Ptr >= 999) M.Hit_Array_Ptr=999;
		if(Tag.Skip) {Gap=0;} else Gap=Tag.End-Tag.Start;
		if (!Gap)
		{
			if (Tag.Skip) Tag.Start = Conversion_Factor-revfmi->saValue[Tag.End/revfmi->saInterval]+Tag.Skip-1;
			else Tag.Start=Conversion_Factor-BWTSaValue(revfmi,Tag.Start);
			Tag.End=Tag.Start;
		}
		M.Hit_Array[M.Hit_Array_Ptr++]=Tag;
		Gap++;M.Hits += Gap;
		M.Stats[In_Mismatch] +=Gap;
	}
}


void Convert_To_Reverse(SARange &Tag,int StringLength,BWT *revfmi,MEMX & M)
{
	char New_Char,Temp_One,Temp_Two;
	char Temp_Char_Array[MAX_MISMATCHES_BOUND];
	char* Current_Tag=M.Current_Tag;
	unsigned pos;
	unsigned Gap=(Tag.Skip)? 0:Tag.End-Tag.Start;
	int c;
	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			pos=Tag.Mismatch_Pos[i];
			Temp_Char_Array[i]=Current_Tag[pos];
			Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
		}

	}
	if(M.Lookupsize == 3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Tag.Start=M.Forward_Start_LookupX[c];Tag.Level=M.Lookupsize+1;
	Tag.Skip=0;
	if(!Tag.Skip)
	{
		while (Tag.Level <= StringLength)
		{
			New_Char=Current_Tag[1-2+Tag.Level];
			Tag.Start = revfmi->cumulativeFreq[New_Char] + BWTOccValue(revfmi, Tag.Start, New_Char) + 1;
			Tag.Level++;
		}
		Tag.End=Tag.Start+Gap;
	}

	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			Current_Tag[Tag.Mismatch_Pos[i]]=Temp_Char_Array[i];
		}
	} 
}
//}-----------------------------------  PRINT ROUTINES ---------------------------------------------------------


//{-----------------------------------  Decode ---------------------------------------------------------
void Open_BAT_File(char* INPUTFILE, FILE* & Data_File,In_File & I)
{
	Header Head;
	I.NORMAL_TAGS=TRUE;
	Data_File=File_Open(INPUTFILE,"rb");
	fseek(Data_File, 0L, SEEK_END);
	I.File_Size = ftello64(Data_File);
	rewind(Data_File);
	//gzread(Data_File,&Head,sizeof(Header));
	if(!fread(&Head,sizeof(Header),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
	if(!(Head.ID[0]=='B'&&Head.ID[1]=='A'&&Head.ID[2]=='T')) {printf("Not a BAT file\n");exit(0);};
	I.MAXHITS=Head.MAXHITS;
	I.STRINGLENGTH=Head.Tag_Length;
	//PRINT_DESC=Head.Print_Desc;
	I.FILETYPE=Head.FILETYPE;
	/*gzread(Data_File,&I.MAX_MISMATCHES,sizeof(I.MAX_MISMATCHES));if (I.MAX_MISMATCHES >5) I.Stat_Size=7*sizeof(unsigned short); else I.Stat_Size=(I.MAX_MISMATCHES+1)*sizeof(unsigned short);
	gzread(Data_File,&I.TAG_COPY_LEN,sizeof(int));
	gzread(Data_File,&I.ROLLOVER,sizeof(char));
	gzread(Data_File,&I.SCANBOTH,sizeof(char));*/
	if(!fread(&I.MAX_MISMATCHES,sizeof(I.MAX_MISMATCHES),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");if (I.MAX_MISMATCHES >5) I.Stat_Size=7*sizeof(unsigned short); else I.Stat_Size=(I.MAX_MISMATCHES+1)*sizeof(unsigned short);
	if(!fread(&I.TAG_COPY_LEN,sizeof(int),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
	if(!fread(&I.ROLLOVER,sizeof(char),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
	if(!fread(&I.SCANBOTH,sizeof(char),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
	if (Head.HITMODE == PAIREND)
	{
		I.NORMAL_TAGS=FALSE;
		if(!fread(&I.Length_Array[1],sizeof(int),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
		if(!fread(&I.Length_Array[2],sizeof(int),1,Data_File)) printf("Open_BAT_File():Error reading input file...\n");
	}
	else {I.Length_Array[1]=I.STRINGLENGTH;I.Length_Array[2]=0;}
}

void Location_To_Genome(unsigned & Location,Ann_Info & A)
{
	std::map <unsigned,Ann_Info> ::iterator I;
	I=Annotations.lower_bound(Location);
	if (I->first != Location) I--;
	A=I->second;	
	//Location=Location-I->first;
	Location=Location-A.Cumulative_Size;
}
/*int  Location_To_Genome(unsigned & Location,unsigned *Offsets,int Genome_Count)
{
	int Genome_Position=0;
	while ( Genome_Position< Genome_Count )
	{
		if (Location < Offsets[Genome_Position]) break;
		Genome_Position++;
	}
	Genome_Position--;
	Location=Location-Offsets[Genome_Position];
	return Genome_Position;
}*/
void Get_Bases2 (unsigned char* Original_Text,unsigned Location,int StringLength,char* Org_String) 
{
	//if (--Location<0) Location=0;
	if(Location >= File_size || (--Location)<0) Location=0;
	//if (unsigned(--Location)<0) Location=0;
	assert (StringLength<SW_STRING_BUFFER);
	{
		for (int i=0;i<=StringLength;i++)
		{
			unsigned char L= (unsigned char) (Original_Text[(Location+i)/4] << (((Location+i) % 4) * 2)) >>6;
			Org_String[i]=L;
		}
	}
}
unsigned Process_GIS(READ & RawR,int & Top_Penalty,int & OrignHits,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Final_Alignments,int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char hitType,std::string readString,char FILETYPE,BWT* fwfmi,BWT* revfmi,Output_Record Record,Mismatches_Record_GIS MismatchesGIS,char* New_Record,int StringLength,char USELOCATION,char PLUSSTRAND,unsigned MAXHITS,unsigned Offset,char* & Output,unsigned* Location_Array,int Genome_Count,Offset_Record* Genome_Offsets,int* Length_Array,char* N,unsigned Hits)
{
	//char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	char* Quality[MAXTAG+1];
	unsigned Conversion_Factor;
	char* Mismatch_Desc_Ptr;
	char* Buffer_Index=Output;
	char Mismatches_Desc[1000];
	unsigned Location;
	int Genome_Position;
	unsigned int minHits=0;
	Ann_Info A;
	Alignment F;
	
	int Mismatch_Count = Record.Mismatches;
	//StringLength=0;
	//for(;readString[StringLength]!=0 && readString[StringLength]!='\n' && readString[StringLength]!='\r';StringLength++);
	StringLength=Length_Array[New_Record[0]];
	Conversion_Factor=revfmi->textLength-StringLength;
	//Ins_Format=(char*)"%d<%c\t";Del_Format=(char*)"%d>D\t";Mis_Format=(char*)"%d>%c\t"; 

	//--------------------------- Decode Mismatches ----------------------------------------------
	Mismatch_Desc_Ptr=Mismatches_Desc;

	*Mismatch_Desc_Ptr=0;
	//--------------------------- Decode Mismatches ----------------------------------------------
	time_t Start_Time,End_Time;//wrong
	time(&Start_Time);//wrong
	for (int j=0;j<=Record.Gap && /*OrignHits*/ Hits < maxhits-1  ;j++) //MAXHITS //&& Hits< maxhits-1 //
	{
		if( Record.Index)//print record...
		{
			if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
			else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
		}
		else
		{
			if (Record.Skip) Location = fwfmi->saValue[Record.Start/fwfmi->saInterval]-Record.Skip+1;
			else Location=BWTSaValue(fwfmi,Record.Start);
		}

		Location -= Offset;
		if (!PLUSSTRAND && New_Record[1]=='-') Location=Location+StringLength;
		Location_To_Genome(Location,A);
		
		char rawRef[400];long Locate=0;
		for(int i=0;i<=String_Hash[A.Name];i++)
		{
			Locate+=Genome_Offsets[i].Offset;
		}
		Locate+=Location;Locate++;
		Get_Bases2(Original_Text_Ori,Locate,StringLength,rawRef);//momo
		
		if(MAX_MISMATCHES>0)
		{
			READ RawRs;
			if(New_Record[1]=='-')
			{
				int i;int Lens=RawR.Real_Len;
				for (i=0;i<=Lens-1;i++){RawRs.Raw_Tag_Copy[Lens-1-i]=Char_To_CharC[RawR.Raw_Tag_Copy[i]];};
			}
			else
			{
				RawRs=RawR;
			}
			char* rawRead =RawRs.Raw_Tag_Copy;

			Mismatch_Count=0;
				
				{
					for(int i=0;i<StringLength;i++)
					{
						assert(rawRef[i]<4 && *rawRef>=0);
						assert(*rawRead!='\0');
						//printf("%u",rawRef[i]);
						if (hitType=='1' || hitType=='3') {
							if (!(*rawRead=='C' && rawRef[i]==1) && !(*rawRead=='T' && rawRef[i]==1) &&  (Char_To_Code[*rawRead] != rawRef[i]) )
							{
								Mismatch_Count++;
							}
						}
						else if (hitType=='2' || hitType=='4') {
							if (!(*rawRead=='G' && rawRef[i]==2) && !(*rawRead=='A' && rawRef[i]==2) && (Char_To_Code[*rawRead] != rawRef[i]) )
							{
								Mismatch_Count++;
							}
						}
						
						rawRead++;
					}
					if(Mismatch_Count<minTrueMis) {Hits=0;minTrueMis=Mismatch_Count;Multipe=false;SecondMisN=0;}
					else if(Mismatch_Count==minTrueMis) Multipe=true;
					
					if(Mismatch_Count>maxTrueMis) {maxTrueMis=Mismatch_Count;}
				}
				/**/
		}
		char CIGr[100];sprintf(CIGr,"%dM",StringLength);
		if(/*Mismatch_Count == minTrueMis  &&*/ Mismatch_Count <= MAX_MISMATCHES) 
		{
			F.source=hitType;F.Chr=A.Name;F.Sign=New_Record[1];F.Loc=Location+1;strcpy(F.Cigar,CIGr);F.QScore=60;F.Mismatch=Mismatch_Count;F.Top_Penalty=Top_Penalty;
			Final_Alignments.push(F);
			//Output+=sprintf(Output,"%c\t%s\t%c\t%u\t%d\t%d\t%s\n",New_Record[0]+'0',A.Name,New_Record[1],Location+1,Mismatch_Count,Length_Array[New_Record[0]],CIGr);//,Record.Tag,Mismatches_Desc
			if( Mismatch_Count == minTrueMis ) Hits++;
			if(minTrueMis==0) minHits++;
		}else if(Mismatch_Count == minTrueMis+1)
		{
			SecondMisN++;
		}
		Record.Start++;OrignHits++;
		if(minTrueMis==0 && minHits >= 10) break;
		time(&End_Time);//wrong
		if(difftime(End_Time,Start_Time) > 2) break;//wrong 120
	}

	return Hits;
//--------------------------- Decode Mismatches ----------------------------------------------
}
/*
unsigned Process_GIS(char FILETYPE,BWT* fwfmi,BWT* revfmi,Output_Record Record,Mismatches_Record_GIS MismatchesGIS,char* New_Record,int StringLength,char USELOCATION,char PLUSSTRAND,unsigned MAXHITS,unsigned Offset,char* & Output,unsigned* Location_Array,int Genome_Count,Offset_Record* Genome_Offsets,int* Length_Array,char* N)
{
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	char* Quality[MAXTAG+1];
	unsigned Conversion_Factor;
	char* Mismatch_Desc_Ptr;
	char* Buffer_Index=Output;
	char Mismatches_Desc[1000];
	unsigned Hits=0;
	unsigned Location;
	int Genome_Position;
	Ann_Info A;
	

	int Mismatch_Count = Record.Mismatches;
	StringLength=Length_Array[New_Record[0]];
	Conversion_Factor=revfmi->textLength-StringLength;
	Ins_Format=(char*)"%d<%c\t";Del_Format=(char*)"%d>D\t";Mis_Format=(char*)"%d>%c\t"; 

	//--------------------------- Decode Mismatches ----------------------------------------------
	Mismatch_Desc_Ptr=Mismatches_Desc;
	if (New_Record[2])
	{
		for(int i=0;i<New_Record[2]*2;i+=2)
		{
			Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,"%d:%c\t",N[i],N[i+1]);
		}
	}

	if (Mismatch_Count)//mismatch positions..
	{
		for( int i=0;i<Mismatch_Count;i++) Mismatch_Desc_Ptr += sprintf(Mismatch_Desc_Ptr,Mis_Format,MismatchesGIS.Mismatch_Pos[i],Code_To_Char[MismatchesGIS.Mismatch_Char>>(2*i) & 3]);
	}
	*Mismatch_Desc_Ptr=0;
	//--------------------------- Decode Mismatches ----------------------------------------------
	for (int j=0;j<=Record.Gap && Hits< MAXHITS;j++)
	{
		if( Record.Index)//print record...
		{
			if (Record.Skip) Location = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
			else Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
		}
		else
		{
			if (Record.Skip) Location = fwfmi->saValue[Record.Start/fwfmi->saInterval]-Record.Skip+1;
			else Location=BWTSaValue(fwfmi,Record.Start);
		}

		Location -= Offset;
		if (!PLUSSTRAND && New_Record[1]=='-') Location=Location+StringLength;
		Location_To_Genome(Location,A);
		if (Location+StringLength > A.Size) 
		{
//Bad_Hit
		}
		Output+=sprintf(Output,"%c\t%s\t%c\t%u\t%d\t%d\t%d\t%s\n",New_Record[0]+'0',Genome_Offsets[Genome_Position].Genome,New_Record[1],Location,Mismatch_Count,Length_Array[New_Record[0]],Record.Tag,Mismatches_Desc);
		Hits++;Record.Start++;
	}

	return Hits;
//--------------------------- Decode Mismatches ----------------------------------------------
}
*/

int Load_Location(char* LOCATIONFILE, std::map <unsigned, Ann_Info> & Annotations,std::map <unsigned, Ann_Info> ::iterator & S,std::map <unsigned, Ann_Info> ::iterator & E,unsigned* Location_Array)
{
	char Genome_Name_Buf[300];
	FILE* Location_File=File_Open(LOCATIONFILE,"r");
	char* Genome_Name;unsigned Off_Cum=0,Off=0;
	int Genome_Count=0;
	while (fgets(Genome_Name_Buf,300,Location_File)!=0)
	{
		Off=atoi(Genome_Name_Buf);
		if (Genome_Count)
		{
			Annotations[Off_Cum].Size=Off;
			Annotations[Off_Cum].Name=Genome_Name;
			Annotations[Off_Cum].Cumulative_Size=Off_Cum;
		}
		Off_Cum+=Off;

		fgets(Genome_Name_Buf,300,Location_File);
		for(int i=0;i<40;i++) 
		{
			if (Genome_Name_Buf[i] == '\n' ||Genome_Name_Buf[i] == '\r')
			{ 
				Genome_Name_Buf[i]=0;
				Genome_Name=new char[i+1];strcpy(Genome_Name,Genome_Name_Buf);
				break;
			} 
		}
		Genome_Count++;	
	}
//-------------
	/*FILE* Location_File=File_Open(LOCATIONFILE,"r");
	int Genome_Count=0;
	while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
	{
		Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
		if(!fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)) break;//printf ("Load_Location():Error reading location file...\n");
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
	for ( int i=1;i<Genome_Count;i++)
	{
		Location_Array[i]=Location_Array[i-1]+Genome_Offsets[i].Offset;
	}
	return Genome_Count;*/
}

int Load_LocationN(char* NLOCATIONFILE, std::map <unsigned, Ann_Info> & Annotations,std::map <unsigned, Ann_Info> ::iterator & S,std::map <unsigned, Ann_Info> ::iterator & E,unsigned* Location_Array)
{
	char Genome_Name_Buf[300];
	FILE* Location_File=File_Open(NLOCATIONFILE,"r");
	char* Genome_Name;unsigned Off_Cum=0,This_Size=0,Size,This_Off=0,Next_Off;
	int Genome_Count=0;
	std::map <unsigned, Ann_Info> TAnnotations=Annotations;
	S=TAnnotations.begin();
	assert(S != TAnnotations.end()); 
	S++;Next_Off=S->first; 

	while (fgets(Genome_Name_Buf,300,Location_File)!=0)
	{
		Size=atoi(Genome_Name_Buf);
		if (Genome_Count)
		{
			This_Size+=Size;
			Annotations[Off_Cum].Size=This_Size;
			Annotations[Off_Cum].Name=Genome_Name;
			Annotations[Off_Cum].Cumulative_Size=This_Off;
		}
		Off_Cum+=Size;
		if (Next_Off && Off_Cum == Next_Off) 
		{
			This_Off=Next_Off;This_Size=0;
			if (++S==TAnnotations.end()) Next_Off=0;else Next_Off=S->first;
		}
		assert(!Next_Off || Off_Cum < Next_Off);

		fgets(Genome_Name_Buf,300,Location_File);
		for(int i=0;i<40;i++) 
		{
			if (Genome_Name_Buf[i] == '\n' ||Genome_Name_Buf[i] == '\r')
			{ 
				Genome_Name_Buf[i]=0;
				Genome_Name=new char[i+1];strcpy(Genome_Name,Genome_Name_Buf);
				break;
			} 
		}
		Genome_Count++;	
	}

	/*
//Print Decode info...
	Ann_Info A;
	S=Annotations.begin();
	for (;S!=Annotations.end();S++)
	{
		A=S->second;
		printf ("%s\n%u:%u:%u\n",A.Name,S->first,A.Cumulative_Size,A.Size); 
	}
	exit(-1);*/
//-------------
	/*FILE* Location_File=File_Open(LOCATIONFILE,"r");
	int Genome_Count=0;
	while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
	{
		Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
		if(!fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)) break;//printf ("Load_Location():Error reading location file...\n");
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
	for ( int i=1;i<Genome_Count;i++)
	{
		Location_Array[i]=Location_Array[i-1]+Genome_Offsets[i].Offset;
	}
	return Genome_Count;*/
}
//}-----------------------------------  Decode ---------------------------------------------------------


//{-----------------------------------  Progress  ---------------------------------------------------------

void Show_Progress(unsigned Percentage)
{
	if (Percentage >=97) return;
	fprintf(stderr,"+%u%%\b\b\b",Percentage);
	fflush(stdout);
}
//}-----------------------------------  Progress  ---------------------------------------------------------

SARange Seed_ExtendF(const char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* & revfmi,MEMX & M)
{
	unsigned Branch_Characters[4];
	if (!Tag.Start) return Tag;
	Start=Start-2;//Adjust for offset difference
	int FSStack_Top=0;
	int Least_Mis=INT_MAX;SARange Least_Hit;bool Unique=false;
	M.FSSStack[0]=Tag;
	SARange Branch_Ranges[4];
	struct SARange Range,Temp_Range;

	while(FSStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.FSSStack[FSStack_Top];
		FSStack_Top--;		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			if(Seed_Forwards_OneSA(Current_Tag,Range,5,Start,StringLength,revfmi,M,Least_Mis,Unique))
			{
				if (Least_Mis >= Range.Mismatches)
				{
					if(Least_Mis==Range.Mismatches) Unique=false;
					else
					{
						Least_Hit=Range;
						Least_Mis=Range.Mismatches; Unique=true;
					}
				}
			}

			//if(MAXHITS<=M.Hits) return;
		}
		else
		{
			Branch_Detect(Range,revfmi,Start,Branch_Characters,Branch_Ranges);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
						Temp_Range.Mismatches++;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches <= Count) //dont print exact matches
							{
								if (Least_Mis >= Temp_Range.Mismatches)
								{
									if(Least_Mis==Temp_Range.Mismatches) Unique=false;
									else
									{
										Least_Hit=Temp_Range;
										Least_Mis=Temp_Range.Mismatches; Unique=true;
									}
								}
								//if (MAXHITS<=M.Hits) return;
							}
							else continue;
						}
						else 
						{

							FSStack_Top++;//Push range
							Temp_Range.Level++;
							M.FSSStack[FSStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	if(!Unique) {Least_Hit.Start=0;}
	return Least_Hit;
}

bool Seed_Forwards_OneSA(const char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M,int & Least_Mis,bool & Unique)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return false;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		//if (!Do_Branch[Tag.Level+Start] && Current_Tag[Tag.Level+Start]!=Now) return;  
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				return true;
			}
			else {Tag.Level++;continue;}
		} 
		else return false;
	}
}

SARange Seed_ExtendB(const char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* & fmi,MEMX & M)
{
	if (!Tag.Start) return Tag;
	unsigned Branch_Characters[4];
	SARange Branch_Ranges[4];
	int BMStack_Top=0;
	M.BMStack[0]=Tag;
	int Least_Mis=M.Least_Mis;SARange Least_Hit;bool Unique=false;
	struct SARange Range,Temp_Range;
	Least_Hit.Start=0;
	while(BMStack_Top!=-1)//While Stack non-empty....
	{
		Range=M.BMStack[BMStack_Top];
		BMStack_Top--;	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			if(Seed_Backwards_OneSA(Current_Tag,Range,Least_Mis,Start,StringLength,fmi,M,Least_Mis,Unique))
			{
				if (Least_Mis >= Range.Mismatches)
				{
					if(Least_Mis==Range.Mismatches) Unique=false;
					else
					{
						Least_Hit=Range;
						Least_Mis=Range.Mismatches; Unique=true;
					}
				}
			}
		}
		else
		{
			Branch_Detect_Backwards(Current_Tag,Range,fmi,Start,Branch_Characters,Branch_Ranges);
			//for(int Branch=3;Branch>=0;Branch--)
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches <= Count) //dont print exact matches
							{
								if (Least_Mis >= Temp_Range.Mismatches)
								{
									if(Least_Mis==Temp_Range.Mismatches) Unique=false;
									else
									{
										Least_Hit=Temp_Range;
										Least_Mis=Temp_Range.Mismatches; Unique=true;
									}
								}
							}
							else continue;
							/*if(Temp_Range.Mismatches==Count)
							{
								Print_LocationX(Temp_Range,M);
								if(MAXHITS<=M.Hits) return;
							}
							else continue;*/
						}
						else
						{
							BMStack_Top++;//Push range
							Temp_Range.Level++;
							M.BMStack[BMStack_Top]=Temp_Range;
						}
					}
				} 
			}
		}
	}
	/*if(!Unique) {Least_Hit.Start=0;}
	else*/ (M.Least_Mis=Least_Mis);
	return Least_Hit;
}

bool Seed_Backwards_OneSA(const char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M,int & Least_Mis,bool & Unique)
{
	unsigned long Index,Now;
	/*char Now_Dbg[MAXDES];int Now_Dbg_Ptr=0;
	char CT_Dbg[MAXDES];int CT_Dbg_Ptr=0;*/

	if (Tag.Start==0) return false;
	if (Count==INT_MAX) Count=LEAST_MIS;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		/*Now_Dbg[Now_Dbg_Ptr++]="ACGT"[Now];
		CT_Dbg[CT_Dbg_Ptr++]="ACGT"[Current_Tag[Start-Tag.Level]];*/
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;
		
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				//Now_Dbg[Now_Dbg_Ptr]=0;
				//CT_Dbg[CT_Dbg_Ptr]=0;
				//printf("%s\n%s\n",Now_Dbg,CT_Dbg);
				return true;
			}
			else {Tag.Level++;continue;}
		} 
		else 
		{
			/*Now_Dbg[Now_Dbg_Ptr]=0;
			CT_Dbg[CT_Dbg_Ptr]=0;
			printf("%s\n%s\n",Now_Dbg,CT_Dbg);*/
			return false;
		}
	}
}


void Convert_To_Fwd(char* Current_Tag,struct SARange & Tag,int StringLength,BWT *fwfmi,MEMX & M)
{	
	unsigned Temp=0;
	int c;
	char New_Char;
	char Mismatch_Count=Tag.Mismatches;
	unsigned pos;
	for( int i=0;i<Mismatch_Count;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Temp=Temp | (Current_Tag[pos]<<i*2);
		Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
	}
	//Temp_BCR[0]=Branch_Characters[0];Temp_BCR[1]=Branch_Characters[1];Temp_BCR[2]=Branch_Characters[2];Temp_BCR[3]=Branch_Characters[3];
	{
		if(M.Lookupsize==3)
		{
			c=Current_Tag[StringLength-1-0] | (Current_Tag[StringLength-1-1]<<2) | (Current_Tag[StringLength-1-2]<<4);// | (Current_Tag[STRINGLENGTH-1-3]<<6) | Current_Tag[STRINGLENGTH-1-4]<<8 | (Current_Tag[STRINGLENGTH-1-5]<<10);//Use lookup table...
		}
		else
		{
			c=Current_Tag[StringLength-1-0] | (Current_Tag[StringLength-1-1]<<2) | (Current_Tag[StringLength-1-2]<<4) | (Current_Tag[StringLength-1-3]<<6) | Current_Tag[StringLength-1-4]<<8 | (Current_Tag[StringLength-1-5]<<10);//Use lookup table...
		}
		Tag.Start=M.Backward_Start_LookupX[c];Tag.End=M.Backward_End_LookupX[c];
		Tag.Level=M.Lookupsize + 1;Tag.Skip=0;
		Search_Backwards_Exact(Current_Tag,Tag,StringLength,StringLength,fwfmi,M);//Backward scan for ?|0
	}
	//Branch_Characters[0]=Temp_BCR[0];Branch_Characters[1]=Temp_BCR[1];Branch_Characters[2]=Temp_BCR[2];Branch_Characters[3]=Temp_BCR[3];
	for( int i=0;i<Tag.Mismatches;i++)
	{
		pos=Tag.Mismatch_Pos[i];
		Current_Tag[pos]=(Temp>>(2*i)) & 3;
	}
	return;
}


int Scan(bool batmeth1,MEMX & MF,MEMX & MC,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int & Top,int Max_Hits)
{
	if(Next_Mis== -1) 
	{
		MF.Hit_Array[MF.Hit_Array_Ptr].Start=0;MC.Hit_Array[MC.Hit_Array_Ptr].Start=0;//tag sentinels to sa lists..
		MF.Hit_Array_Ptr++;MC.Hit_Array_Ptr++;//Setup for suboptimal hits..
		return -1;
	}
	assert(Next_Mis >=0);assert(MAX_MISMATCHES >= Next_Mis);assert (Next_Mis <= 5);
	assert(MAX_MISMATCHES<=Max_MM);
	//if(BWAPOL && Next_Mis) MAX_MISMATCHES = (Next_Mis<MAX_MISMATCHES) ? Next_Mis : MAX_MISMATCHES;
	if(Next_Mis) MAX_MISMATCHES = (Next_Mis<MAX_MISMATCHES) ? Next_Mis : MAX_MISMATCHES;
	int In_Mis=0,Hits=0;MF.Hits=MC.Hits=0;
	if (Next_Mis == 0) goto Zero; else if (Next_Mis ==1) goto One;else if (Next_Mis ==2) goto Two;else if (Next_Mis ==3) goto Three;else if (Next_Mis ==4) goto Four; else goto Five;
Zero:
	Hits+=Zero_Mismatch(batmeth1,MF.Current_Tag,L,revfmi,MF,fwfmi);
	Hits+=Zero_Mismatch(batmeth1,MC.Current_Tag,L,revfmi,MC,fwfmi);
One:
	if (!Hits && MAX_MISMATCHES >0)
	{
		In_Mis=1;
		Hits+=One_Mismatch(batmeth1,MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		Hits+=One_Mismatch(batmeth1,MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Two:
	if (!Hits && MAX_MISMATCHES >1)
	{
		In_Mis=2;
		Hits+=Two_Mismatch(batmeth1,MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		Hits+=Two_Mismatch(batmeth1,MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Three:
	if (!Hits && MAX_MISMATCHES >2)
	{
		In_Mis=3;
		Hits+=Three_Mismatch(batmeth1,MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		Hits+=Three_Mismatch(batmeth1,MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Four:
	if (!Hits && MAX_MISMATCHES >3)
	{
		In_Mis=4;
		Hits+=Four_Mismatch(batmeth1,MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		Hits+=Four_Mismatch(batmeth1,MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Five:
	if (!Hits && MAX_MISMATCHES >4)
	{
		In_Mis=5;
		Hits+=Five_Mismatch(batmeth1,MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		Hits+=Five_Mismatch(batmeth1,MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}

	MF.Hit_Array[MF.Hit_Array_Ptr].Start=0;MC.Hit_Array[MC.Hit_Array_Ptr].Start=0;//tag sentinels to sa lists..
	MF.Hit_Array_Ptr++;MC.Hit_Array_Ptr++;//Setup for suboptimal hits..
	assert(In_Mis <= MAX_MISMATCHES);
	Top=Hits;
	return (Hits ? In_Mis : -1) ;
}
//----------------------------batmeth--

void Print_Hits(READ & RawR,int & Last_Mis,bool & batmeth1,std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> & Final_Alignments,int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char type,Offset_Record *Genome_Offsets,std::string readString,MEMX & M,LEN & L,FILE* Output_File, char ONEFMINDEX,BWT *revfmi,BWT *fwfmi, OUTPUT & O)
{
	Output_Record Record;
	char New_Record[4];
	char* Output=O.Buffer;
	Mismatches_Record_GIS MismatchesGIS;

	//>>Output+=sprintf(Output,"@\n");
	int Desc_End;for(Desc_End=0;M.Read.Description[Desc_End]!='\n' && M.Read.Description[Desc_End]!='\r' && M.Read.Description[Desc_End]!=0;Desc_End++);M.Read.Description[Desc_End]=0;
	M.Read.Tag_Copy[L.STRINGLENGTH]=0;M.Read.Quality[L.STRINGLENGTH]=0;
	//>>Output+=sprintf(Output,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\n",M.Read.Description,M.Read.Tag_Copy,M.Read.Quality,M.Stats[0],M.Stats[1],M.Stats[2],M.Stats[3],M.Stats[4],M.Stats[5],M.Stats[6]);
	int Hits=0;
	int Top_Penalty=0;
	int i=0;int OrignHits=0;
	while (i <M.Hit_Array_Ptr)
	{
		SARange Tag=M.Hit_Array[i++];

		if(Last_Mis==0)
		{
			SARange SA=M.Exact_Match_Forward[L.LH-1];
			if(SA.Start)//Get stat for later mapQ calculation..
			{
				if (SA.Skip)
					Top_Penalty++;
				else
				Top_Penalty=SA.End-SA.Start;
			}
		}
	
		if ( ONEFMINDEX && FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse_batmeth(batmeth1,Tag,L,revfmi,M,fwfmi);
			if (Tag.Skip) Tag.Start=Tag.End;
		}


		unsigned Gap;
		if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
		Record.Gap=Gap;


		//New_Record[0]='%';
		New_Record[0]=M.Read.Tag_Number;
		New_Record[1]=Tag.Strand;
		New_Record[2]=M.Read.NCount;
		//fwriteX(&New_Record,4,1,Output_File);

		memcpy(MismatchesGIS.Mismatch_Pos,Tag.Mismatch_Pos,Tag.Mismatches);//MAX_MISMATCHES_BOUND);
		MismatchesGIS.Mismatch_Char=Tag.Mismatch_Char;
		Record.Start=Tag.Start;
		Record.Tag=M.Read.Read_Number;
		Record.Skip=Tag.Skip;
		Record.Mismatches=Tag.Mismatches;
		if (ONEFMINDEX) Record.Index=REVERSE; else Record.Index=Tag.FMIndex;
		int Skip_Length=Tag.Mismatches+sizeof(unsigned);
		//fwriteX(&Record,sizeof(Record),1,Output_File);
		//fwriteX(&MismatchesGIS,Skip_Length,1,Output_File);//;Write_Buffer_Ptr += Skip_Length;//+1;
		Hits=Process_GIS(RawR,Top_Penalty,OrignHits,Final_Alignments,SecondMisN,Multipe,minTrueMis,maxTrueMis,type,readString,O.FILETYPE,fwfmi,revfmi,Record,MismatchesGIS,New_Record,L.STRINGLENGTH,TRUE,O.PLUSSTRAND,O.MaxHits,O.Offset,Output,O.Location_Array,O.Genome_Count,Genome_Offsets,O.Length_Array,M.Read.N,Hits);
	}
	O.Buffer_End=Output;
}

void Convert_To_Reverse_batmeth(bool batmeth1,SARange &Tag,LEN & L,BWT *revfmi,MEMX & M,BWT *fwfmi)
{
	char New_Char,Temp_One,Temp_Two;
	char Temp_Char_Array[MAX_MISMATCHES_BOUND];
	char* Current_Tag=M.Current_Tag;
	unsigned pos;
	unsigned Gap=Tag.End-Tag.Start-1;
	int c;

	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			pos=Tag.Mismatch_Pos[i];
			Temp_Char_Array[i]=Current_Tag[pos];
			Current_Tag[pos]=Tag.Mismatch_Char>>(2*i) & 3;
		}

	}
	if(M.Lookupsize == 3)
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4);// | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	else
	{
		c=Current_Tag[0] | (Current_Tag[1]<<2) | (Current_Tag[2]<<4) | (Current_Tag[3]<<6) | Current_Tag[4]<<8 | (Current_Tag[5]<<10);//Use lookup table...
	}
	Tag.Start=M.Forward_Start_LookupX[c];Tag.Level=M.Lookupsize+1;
	Tag.Skip=0;
	if(!Tag.Skip)
	{
		while (Tag.Level <= L.STRINGLENGTH)
		{
			New_Char=Current_Tag[1-2+Tag.Level];
			Tag.Start = revfmi->cumulativeFreq[New_Char] + BWTOccValue(revfmi, Tag.Start, New_Char) + 1;
			Tag.Level++;
		}
		Tag.End=Tag.Start+Gap;
	}
	else
	{
		Search_Forwards_Exact(batmeth1,Current_Tag,Tag,-1,L.STRINGLENGTH,revfmi,M,L.LH,revfmi,fwfmi);
		/*Tag.Skip=0;
		while (Tag.Level <= STRINGLENGTH)
		{
			Get_SARange_Fast(Current_Tag[Tag.Level-1],Tag,revfmi);
			if (Tag.Skip) Tag.Skip++;
			else if(Tag.End==Tag.Start && Tag.Start % SAINTERVAL == 0) 
			{
				Tag.Skip++;Tag.End=Tag.Start;
			}
		}*/
	}

	if(Tag.Mismatches)
	{
		for( int i=0;i<Tag.Mismatches;i++)
		{
			Current_Tag[Tag.Mismatch_Pos[i]]=Temp_Char_Array[i];
		}
	} 
}
