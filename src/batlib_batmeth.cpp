//using namespace std;
#include <stdio.h>
#include <cassert>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <signal.h> 
//#include <dvec.h>
#include <getopt.h>
#include "zlib.h"
#include <queue>
#include <map>
#include <ctype.h>

#include "batlib_batmeth.h"
extern unsigned Genome_Count;
extern std::map <unsigned, Ann_Info> Annotations;
extern std::map <std::string,int> String_Hash;
extern Gene_Hash Genome_List[80];
extern int MAX_MISMATCHES;
extern int minMis;
extern int maxhits;
extern int allHits;
extern int indelMis;

char OUTPUT_ZIPPED;
char SOLID;
char NPOLICY;
char Random_Array[]="tatacgataggacaatgtcttcgaagcccacgcggtaagccggtcattgcggttgtgcgaacactatcagcctcgctgcatggttaccctgggtggataggacgtttgcccgacattttgacacgcataaaaggtctgtagtgggggtggcacaccataaaccctggggcggctccacgatcgtaaaatcctgcgatctg";
const char* Code_To_Char="acgt";const char* Code_To_CharCAPS="ACGT";
char Char_To_CodeC[255];
char Char_To_Code[255];
//{-----------------------------------  INIT ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Analyze_File
 *  Description:  Analyse a file and try to determine its type..
 * =====================================================================================
 */
//void Analyze_File(char & FILETYPE,char & SOLID,FILE* Input_File,int & TAG_COPY_LEN,char & TAB_SEPERATED,int &  PAIR_LENGTH_LEFT,int & PAIR_LENGTH_RIGHT,LEN & L)
void GET_LEN(LEN & L,READ & Read)
{
	L.STRINGLENGTH=strlen(Read.Tag_Copy)-1;
	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
}
void Initial_Length(LEN & L)
{
	L.STRINGLENGTH=70;
	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
}
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
		if(!fgets(Description,MAXDES,I.Input_File)) {printf("Analyze_File(): error reading file...\n");exit(-1);}
		if (Description[0] != '#') break;
		Last=ftello64(I.Input_File);//mark last comment...
	}
	if(!fgets(Current_Tag,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}
	if(Current_Tag[2]>='0' && Current_Tag[2] <='3') I.SOLID=TRUE;else I.SOLID =FALSE;//csfasta has numbers..
	for(I.TAG_COPY_LEN=0;Current_Tag[I.TAG_COPY_LEN]!='\n' && Current_Tag[I.TAG_COPY_LEN]!='\r' && Current_Tag[I.TAG_COPY_LEN]!=0;I.TAG_COPY_LEN++);// Find the length of tag line.should be the same length is quality is present.. 
	for(L.STRINGLENGTH=0;Current_Tag[L.STRINGLENGTH]!='\n' && Current_Tag[L.STRINGLENGTH]!='\r' && Current_Tag[L.STRINGLENGTH]!=0 && Current_Tag[L.STRINGLENGTH]!=PAIR_END_SEPERATOR;L.STRINGLENGTH++);//scan for a split that indicates tab seperated PET
	if(Current_Tag[L.STRINGLENGTH]==PAIR_END_SEPERATOR) 
	{
		I.TAB_SEPERATED=TRUE;//we have pair ended tags..
		I.PAIR_LENGTH_LEFT=L.STRINGLENGTH;
		for(I.PAIR_LENGTH_RIGHT=0;Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\n' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\r' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!=0;I.PAIR_LENGTH_RIGHT++);
	}
	else I.TAB_SEPERATED=FALSE;
	if(!fgets(Quality,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}//plus
	if (Quality[0]=='>') I.FILETYPE=FA;
	else 
	{
		I.FILETYPE=FQ;
		if(Quality[0] != '+' || Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
	}
	fseek(I.Input_File,Last,SEEK_SET);//go top
	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
}

void Analyze_File(INFILE & I)
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
		if(!fgets(Description,MAXDES,I.Input_File)) {printf("Analyze_File(): error reading file...\n");exit(-1);}
		if (Description[0] != '#') break;
		Last=ftello64(I.Input_File);//mark last comment...
	}
	if(!fgets(Current_Tag,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}
	if(Current_Tag[2]>='0' && Current_Tag[2] <='3') I.SOLID=TRUE;else I.SOLID =FALSE;//csfasta has numbers..
	for(I.TAG_COPY_LEN=0;Current_Tag[I.TAG_COPY_LEN]!='\n' && Current_Tag[I.TAG_COPY_LEN]!='\r' && Current_Tag[I.TAG_COPY_LEN]!=0;I.TAG_COPY_LEN++);// Find the length of tag line.should be the same length is quality is present.. 
/*	for(L.STRINGLENGTH=0;Current_Tag[L.STRINGLENGTH]!='\n' && Current_Tag[L.STRINGLENGTH]!='\r' && Current_Tag[L.STRINGLENGTH]!=0 && Current_Tag[L.STRINGLENGTH]!=PAIR_END_SEPERATOR;L.STRINGLENGTH++);//scan for a split that indicates tab seperated PET
	if(Current_Tag[L.STRINGLENGTH]==PAIR_END_SEPERATOR) 
	{
		I.TAB_SEPERATED=TRUE;//we have pair ended tags..
		I.PAIR_LENGTH_LEFT=L.STRINGLENGTH;
		for(I.PAIR_LENGTH_RIGHT=0;Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\n' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!='\r' && Current_Tag[L.STRINGLENGTH+1+I.PAIR_LENGTH_RIGHT]!=0;I.PAIR_LENGTH_RIGHT++);
	}
	else I.TAB_SEPERATED=FALSE;
*/
	if(!fgets(Quality,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}//plus
	if (Quality[0]=='>') I.FILETYPE=FA;
	else 
	{
		I.FILETYPE=FQ;
		if(Quality[0] != '+' || Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
	}
	fseek(I.Input_File,Last,SEEK_SET);//go top
//	L.STRINGLENGTH_ORG=L.STRINGLENGTH;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Lookup_Size
 *  Description:  calulate lookupsize for given parameters...
 * =====================================================================================
 */
int Get_Lookup_Size(char MAX_MISMATCHES,char STRINGLENGTH)
{
	int LOOKUPSIZE = 3;//6;
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
	int Max_Limit=5;//5;//3
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
		printf("Load_Indexes():FM index load error \n"); 
		exit(1);
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
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....\n",File_Name);
		exit(1);
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
		//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
		if(!fgets(Read.Tag_Copy,MAXDES,Input_File)) {printf ("Read_Tag():Error reading file..\n");exit(0);};//tag
		if (FILETYPE == FQ)
		{
			//gzgets(Input_File,Plus,MAXTAG);//plus
			if(!fgets(Read.Plus,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file..\n");exit(0);};//plus
			//gzgets(Input_File,Quality,MAXTAG);//phred
			if(!fgets(Read.Quality,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file..\n");exit(0);};//phred
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
bool Process_Read(READ & R, BATREAD & B,MEMX & MF, MEMX & MC)
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
	if(R.NCount>2) {/*if(PRINT_MISHITS) fprintf(Mishit_File,"%s%s%s%s", Description,TAg_Copy,Plus,Quality);*/ return false;} //jqjq_Hardcode
	//if (NISMISMATCH && MAX_MISMATCHES <= R.NCount ) continue;//Too many N's than mismatches...
	B.Forward[B.StringLength]='+';
	B.Complement[B.StringLength]='-';


	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Read=R;MF.Current_Tag=B.Forward;
	MC.Hits=0;MC.Hit_Array_Ptr=0;MC.Read=R;MC.Current_Tag=B.Complement;
	for (int i=0;i<7;i++) {MF.Stats[i]=0;MC.Stats[i]=0;}
	return true;

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
		if (LRange.Start>LRange.End) LRange.Start=0;
		LRange.Label=i;
		for(int j=0;j<4;j++) Build_Preindex_Forward(revfmi,LRange, 2, j,M);
		LRange.Start=fwfmi->cumulativeFreq[i] + 1;
		LRange.End=fwfmi->cumulativeFreq[i + 1];
		if (LRange.Start>LRange.End) LRange.Start=0;
		LRange.Label=i;
		for(int j=0;j<4;j++) Build_Preindex_Backward(fwfmi,LRange, 2, j,M);

	}
}

void Build_Preindex_Backward_Buggy(BWT* fwfmi,Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		M.Backward_Start_LookupX[Range.Label] = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
		M.Backward_End_LookupX[Range.Label] = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
	}
	else
	{


		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		Range.Start = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
		Range.End = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
		Level ++;
		for ( int i=0;i<4;i++)
		{
			Build_Preindex_Backward( fwfmi,Range, Level,i,M);
		}

	}

}


void Build_Preindex_Backward(BWT* fwfmi,Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if (Range.Start)
		{
			Range.Start = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
			Range.End = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		M.Backward_Start_LookupX[Range.Label] = Range.Start;
		M.Backward_End_LookupX[Range.Label] = Range.End;
	}
	else
	{


		assert (Range.Start<=Range.End);
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
			Build_Preindex_Backward( fwfmi,Range, Level,i,M);
		}

	}

}

void Build_Preindex_Forward(BWT *revfmi, Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		if(Range.Start)
		{
			Range.Start = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
			Range.End = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
			if (Range.Start>Range.End) Range.Start=0;
		}
		M.Forward_Start_LookupX[Range.Label] = Range.Start;
		M.Forward_End_LookupX[Range.Label] = Range.End;
	}
	else
	{


		assert (Range.Start<=Range.End);
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		if(Range.Start)
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

void Build_Preindex_Forward_Buggy(BWT *revfmi, Range Range, int Level, int Bit, const MEMLOOK & M)
{

	if (M.Lookupsize==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		M.Forward_Start_LookupX[Range.Label] = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
		M.Forward_End_LookupX[Range.Label] = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
	}
	else
	{


		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		Range.Start = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
		Range.End = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
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

unsigned Zero_Mismatch(char* Current_Tag,LEN L, BWT* revfmi,MEMX & M)
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
	Search_Forwards_Exact(Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,M,L.LH);//Find exact matches and report... if not found get the range for 0|?
	return M.Hits;

}


//{------------------------------------------- ONE MISMATCH ---------------------------------------------------------------------------------------------
unsigned One_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
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
		Search_Forwards(Current_Tag,Range,1,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for one mismatches of the form 0|1, store possible two mismatches of the form 0|2...
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
		Search_Backwards(Current_Tag,Range,1,L.LH,L.LH,MAXHITS,fwfmi,M);//Backward scan for one mismatches of the form 1|0, store possible mismatches of the form 2|0
		if(MAXHITS<=M.Hits) return M.Hits;
	}
	return M.Hits;

}

void Search_Backwards(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,MEMX & M)
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
			Search_Backwards_OneSA(Current_Tag,Range,Count,Start,StringLength,fmi,M);
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
								Print_LocationX(Temp_Range,M);
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

void Search_Backwards_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M)
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
					Print_LocationX(Tag,M);
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

unsigned Two_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	int c;
	SARange Range,TRange;
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			(M.Two_Mismatches_At_End_Forward[i]).Mismatch_Pos[1]= (L.STRINGLENGTH-1);//mismatches of the form 0|2, with last mismatch at the end...
			Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);
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
			Print_LocationX(M.Two_Mismatches_At_End[i],M);//Mismatches of the form 2|0, with one mismatch at the first position
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
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],2,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for possible two mismatches of the form 0|2, and store candidates for 0|3
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
			Search_Backwards(Current_Tag,M.Mismatches_Backward[i],2,L.LH,L.LH,MAXHITS,fwfmi,M);//scan for possible two mismatches of the form 2|0, and stores the candidates for 3|0
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
	Search_Backwards_X10(Current_Tag,Range,1,L.LH + L.RHQL, L.RHQL,MAXHITS,L,fwfmi,M);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
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
	Search_X01(Current_Tag,TRange,1,L.LH + L.RHQL +1,L.RHQR,MAXHITS,L,fwfmi,revfmi,M);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	return M.Hits;
}

void Search_X01(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT* fwfmi,BWT *revfmi,MEMX & M)
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
			Search_X01_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,L,fwfmi,revfmi,M);
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
									Search_Backwards(Current_Tag,Temp_Range,2,L.LH,L.LH,MAXHITS,fwfmi,M);
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

void Search_X01_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fwfmi,BWT *revfmi,MEMX & M)
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
						Search_Backwards(Current_Tag,Tag,2,L.LH,L.LH,MAXHITS,fwfmi,M);
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

void Search_Backwards_X10(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fwfmi,MEMX & M)
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
			Search_Backwards_X10_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,L,fwfmi,M);
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
								Search_Backwards(Current_Tag,Temp_Range,2,L.LH,L.LH,MAXHITS,fwfmi,M);
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


void Search_Backwards_X10_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fmi,MEMX & M)
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
					Search_Backwards(Current_Tag,Tag,2,L.LH,L.LH,MAXHITS,fmi,M);
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
unsigned Three_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);//mismatches of the form 0|3, with last mismatch at the end...
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
			Print_LocationX(M.Two_Mismatches_At_End[i],M);//Mismatches of the form 3|0, with one mismatch at the first position
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
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],3,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for possible three mismatches of the form 0|3, and finds mismatches of the form 1|2, stores possibles in the form 1|3
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
			Search_Backwards(Current_Tag,M.Mismatches_Backward[i],3,L.LH,L.LH,MAXHITS,fwfmi,M);//scan for possible mismatches of the form 3|0, 2|1 and sotres the candidates for 4|0, 3|1
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
	Search_01X(Current_Tag,Range,1,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M);
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
	Search_10X(Current_Tag,Range,1,L.LHQL, L.LHQL,L,MAXHITS,fwfmi,revfmi,M);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	return M.Hits;
}

void Search_10X_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M)
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
						Search_Forwards(Current_Tag,Tag,3,L.LH+1,L.RH,MAXHITS,revfmi,M);
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
void Search_10X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L, int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M)
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
			Search_10X_OneSA(Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,fmi,revfmi,M);
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
									Search_Forwards(Current_Tag,Temp_Range,3,L.LH+1,L.RH,MAXHITS,revfmi,M);
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
void Search_01X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M)
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
			Search_01X_OneSA(Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,revfmi,M);
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
								Search_Forwards(Current_Tag,Temp_Range,3,L.LH+1,L.RH,MAXHITS,revfmi,M);
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


void Search_01X_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M)
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
					Search_Forwards(Current_Tag,Tag,3,L.LH+1,L.RH,MAXHITS,fmi,M);
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
unsigned Four_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{
	M.FMIndex=REVERSE;
	SARange Range,TRange;

	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);//mismatches of the form 0|4, with last mismatch at the end...
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
			Print_LocationX(M.Two_Mismatches_At_End[i],M);//mismatches of the form 0|4, with one mismatch at the start...
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
		if (M.Possible_04_Pointer > 46000) {M.Mismatches_Forward_Pointer=0;Wrap=TRUE;}
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
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],4,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for possible four mismatches of the form 0|4, and finds mismatches of the form 1|3, stores possibles in the form 1|4
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
			Search_Backwards(Current_Tag,M.Mismatches_Backward[i],4,L.LH,L.LH,MAXHITS,fwfmi,M);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
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
				Search_Forwards(Current_Tag,M.Left_Mishits[i],4,1,L.STRINGLENGTH,MAXHITS,revfmi,M);//find mismatches of the form 022 form, stores possibles of the form 023
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
					Search_Forwards(Current_Tag,M.Right_Mishits[i],4,L.LH+1,L.RH,MAXHITS,revfmi,M);
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
	Search_01LX(Current_Tag,Range,1,LHQLrx +1,LHQRrx,MAXHITS,revfmi,L,M);
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
	Search_10LX(Current_Tag,Range,1,LHQLrx, LHQLrx,MAXHITS,fwfmi,revfmi,L,M);//search for three mismatches of the form 1|2 and stores the candidates for 1|3
	return M.Hits;
}

void Search_10LX_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
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
						Search_Half_Tag_11X(Current_Tag,Tag,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M);
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

void Search_10LX(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
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
			Search_10LX_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,fwfmi,revfmi,L,M);
			if (M.Last_Mismatch_Written ==8 ) return;
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
									Search_Half_Tag_11X(Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M);
									if (M.Last_Mismatch_Written ==8 ) return;
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
void Search_01LX_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,LEN & L,MEMX & M)
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
					Search_Half_Tag_11X(Current_Tag,Tag,2,L.LHQL +1,L.LHQR,L,MAXHITS,fmi,M);
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

void Search_01LX(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *revfmi,LEN & L,MEMX & M)
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
			Search_01LX_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,revfmi,L,M);
			if (M.Last_Mismatch_Written ==8 ) return;
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
								Search_Half_Tag_11X(Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M);
								if (M.Last_Mismatch_Written ==8 ) return;
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


void Search_11X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M)
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
						Search_Half_Tag_11X(Current_Tag,Temp_Range,2,L.LHQL +1,L.LHQR,L,MAXHITS,revfmi,M);
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

void Search_Half_Tag_11X_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M)
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
				Search_Forwards(Current_Tag,Tag,4,L.LH+1,L.RH,MAXHITS,fmi,M);
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

void Search_Half_Tag_11X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M)
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
			Search_Half_Tag_11X_OneSA(Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,revfmi,M);
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
							Search_Forwards(Current_Tag,Temp_Range,4,L.LH+1,L.RH,MAXHITS,revfmi,M);
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
unsigned Five_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M)
{

	SARange Range,TRange;
	M.FMIndex=REVERSE;
	if(M.Two_Mismatches_At_End_Forward_Pointer)//give priority to forward direction as most erros occur in the end..
	{
		for(int i=0;i<M.Two_Mismatches_At_End_Forward_Pointer;i++)
		{
			Print_LocationX(M.Two_Mismatches_At_End_Forward[i],M);//mismatches of the form 0|4, with last mismatch at the end...
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
			Print_LocationX(M.Two_Mismatches_At_End[i],M);//mismatches of the form 0|5, with one mismatch at the start...
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
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],5,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for possible four mismatches of the form 0|5
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;

	}


	M.FMIndex=REVERSE;
	if(M.Mismatches_Forward_Pointer)
	{
		for(int i=M.Mismatches_Forward_Pointer_Last4;i<M.Mismatches_Forward_Pointer_Last5;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],5,1,L.STRINGLENGTH,MAXHITS,revfmi,M);//scan for possible five mismatches of the form 0|5, and finds mismatches of the form 1|4,2|3 
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=REVERSE;
	if(M.Mismatches_Forward_Pointer)
	{
		for(int i=M.Mismatches_Forward_Pointer_Last5;i<M.Mismatches_Forward_Pointer;i++)//Mismatches_Forward_Pointer;i++)
		{
			Search_Forwards(Current_Tag,M.Mismatches_Forward[i],5,L.LH+1,L.RH,MAXHITS,revfmi,M);//scan for possible four mismatches of the form 0|5
			if(MAXHITS<=M.Hits) break;
		}
		if(MAXHITS<=M.Hits) return M.Hits;
	}

	M.FMIndex=FORWARD;
	if(M.Mismatches_Backward_Pointer!=M.Possible_50_Pointer)
	{
		for(int i=M.Mismatches_Backward_Pointer-1;i>=M.Possible_50_Pointer;i--)
		{
			Search_Backwards(Current_Tag,M.Mismatches_Backward[i],5,L.LH,L.LH,MAXHITS,fwfmi,M);//scan for possible mismatches of the form 4|0, 3|1 and sotres the candidates for 5|0, 4|1
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
				Search_Backwards(Current_Tag,M.Possible_20[i],5,L.LH,L.LH,MAXHITS,fwfmi,M);
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
					Search_Backwards(Current_Tag,M.Possible_02[i],5,L.LH,L.LH,MAXHITS,fwfmi,M);
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

	Search_Backwards_XL10(Current_Tag,Range,1,L.LH + L.RHQL+RHQLlx, RHQLlx,MAXHITS,fwfmi,L,M);//?|1|0 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
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
	Search_XL01(Current_Tag,TRange,1,L.LH + L.RHQL+RHQLlx +1,RHQRlx,MAXHITS,fwfmi,revfmi,L,M);//?|0|1 and extend, finds mismatches of the form 1|1 and stres candidates for 2|1
	return M.Hits;
}

void Search_XL01_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* fwfmi,BWT *revfmi,LEN & L,MEMX & M)
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
						Search_Half_Tag_X11(Current_Tag,Tag,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fwfmi,M);
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
void Search_XL01(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M)
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
			Search_XL01_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,fwfmi,revfmi,L,M);
			if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
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
									Search_Half_Tag_X11(Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fwfmi,M);
									if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
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
void Search_Backwards_XL10_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,LEN & L, MEMX & M)
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
					Search_Half_Tag_X11(Current_Tag,Tag,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M);
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
void Search_Backwards_XL10(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *fmi,LEN & L,MEMX & M)
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
			Search_Backwards_XL10_OneSA(Current_Tag,Range,Count,Start,StringLength,MAXHITS,fmi,L,M);
			if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
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
								Search_Half_Tag_X11(Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M);
								if (!M.Larger_Than_Ten && M.Last_Mismatch_Written>5 && M.Last_Mismatch_Written <=10 ) return;
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

void Search_X11(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,unsigned MAXHITS, LEN & L,BWT *fmi,MEMX & M)
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
							Search_Half_Tag_X11(Current_Tag,Temp_Range,2,L.STRINGLENGTH,L.RH,L,MAXHITS,fmi,M);
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

void Search_Half_Tag_X11_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M)
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
				Search_Backwards(Current_Tag,Tag,5,L.LH,L.LH,MAXHITS,fmi,M);//LH,fwfmi);
				return;
			}
			else {Tag.Level++;continue;}
		} 
		else return;
	}
}

void Search_Half_Tag_X11(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M)
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
			Search_Half_Tag_X11_OneSA(Current_Tag,Range,Count,Start,StringLength,L,MAXHITS,fmi,M);
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
							Search_Backwards(Current_Tag,Temp_Range,5,L.LH,L.LH,MAXHITS,fmi,M);//LH,fwfmi);
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
unsigned Guess_Orientation(BWT* fwfmi,BWT* revfmi,MEMX & MF,MEMX & MC,LEN & L,GUESS & G,BATREAD & B)
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
	Search_Forwards_Exact(Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,MF,L.LH);//Find exact matches and report... if not found get the range for 0|?
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
	Search_Forwards_Exact(Current_Tag,Range,-1,L.STRINGLENGTH,revfmi,MC,L.LH);//Find exact matches and report... if not found get the range for 0|?
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
void Search_Forwards_Exact(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi,MEMX & M,int LH)
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
						Print_LocationX(Tag,M);
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
					Print_LocationX(Tag,M);
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

void Search_Forwards(const char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* & revfmi,MEMX & M)
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
			Search_Forwards_OneSA(Current_Tag,Range,Count,Start,StringLength,revfmi,M);
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
								Print_LocationX(Temp_Range,M);
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

void Search_Forwards_OneSA(const char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M)
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
					Print_LocationX(Tag,M);
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

void Print_Hits_To_BAT(MEMX & M,LEN & L,FILE* Output_File, char ONEFMINDEX,BWT *revfmi,char FILETYPE)
{
	Output_Record Record;
	char New_Record[4];
	Mismatches_Record_GIS MismatchesGIS;

	flockfile(Output_File);
	fprintfX(Output_File,(char*)"@%s",M.Read.Description);
	fwriteX(M.Stats,sizeof(unsigned short),M.Stat_Size,Output_File);
	/*if(!NORMAL_TAGS)
	{
		fwriteX(Stats2,sizeof(unsigned short),Stat_Size,Output_File);
	}*/
	fwriteX(M.Read.Tag_Copy,1,L.STRINGLENGTH_ORG,Output_File);//Write_Buffer_Ptr += TAG_COPY_LEN;//write tag...
	if (FILETYPE==FQ) 
	{
		fwriteX(M.Read.Quality,1,L.STRINGLENGTH_ORG,Output_File);
		//if(NORMAL_TAGS) fwriteX(Quality,1,TAG_COPY_LEN,Output_File);
		//else fwriteX(Quality_Copy,1,TAG_COPY_LEN,Output_File);
	}

	int i=0;
	while (i <M.Hit_Array_Ptr)
	{
		SARange Tag=M.Hit_Array[i++];


		if ( ONEFMINDEX && FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse(Tag,L,revfmi,M);
			if (Tag.Skip) Tag.Start=Tag.End;
		}

		/*if (NCount && Tag.Mismatches)
		  {
		//search for places where N's and mismatches overlap
		int Mis_InN=0,Mis_NotinN=0;
		for (int i=0;i<Tag.Mismatches;i++) 
		if (NLocations[Tag.Mismatch_Pos[i]])
		{Mismatches_InN[Mis_InN++]=Current_Tag[Tag.Mismatch_Pos[i]];}
		//If mismatches fall on an N, discard it...
		Mis_NotinN=Tag.Mismatches-Mis_InN;
		if(NISMISMATCH){if (NCount+Mis_NotinN > MAX_MISMATCHES ) return;}
		else if (Mis_NotinN > MAX_MISMATCHES) return;
		}*/

		unsigned Gap;
		if(Tag.Skip) {Gap=0;Tag.Start=Tag.End;} else Gap=Tag.End-Tag.Start;
		Record.Gap=Gap;

		/*if(STATMODE)
		  {
		  if (Stats[0] +Stats[1] > THRESHOLD01 ||Stats[0]+Stats[1]+Stats[2]>THRESHOLD012 )//bad tag....
		  {
		  Write_Buffer_Ptr=Write_Buffer;
		  Total_Tags -= Tags_From_Head;Total_Hits -=Hits1;
		  Hits=MAXHITS;
		  Rollover_Step=TRUE;
		  Tag_Stat_Bad=TRUE;
		  return;
		  }
		  }

		  if ( ONEFMINDEX && FORWARD==FMIndex )// forward search index...
		  {
		  Convert_To_Reverse(Tag);
		  if (Tag.Skip) Tag.Start=Tag.End;
		  }

		  if (!Hits) // First Hit, so write header...
		  {
		  if(!Tag_Printed) Tag_Printed=TRUE;//A hit for this tag has been found...
		  if (SCANBOTH && First_Pass_Hits) {Total_Tags--;Tags_From_Head--;}
		  Total_Tags++;Tags_From_Head++;
		  }

		  if (COUNT_ALLHITS)
		  {
		  if(Hits + Gap > MAXHITS ) Hits = MAXHITS; else Hits=Hits+Gap;//COUNT_ALLHITS =1 => count all the hits, COUNT_ALLHITS =0 => count saranges... i.e. unique hits...
		  }


		  if (HIGHSCAN && Last_Mismatch_Written)//choose best hits...
		  {
		  if (Tag.Mismatches >= Last_Mismatch_Written)//we already have a better hit..
		  {
		  return;
		  }
		  else //we have a better hit... 
		  {
		  Write_Buffer_Ptr=Last_Write_Buffer;
		  }
		  }

		  if ((FILTERUNIQUEHITS && Hits !=1))
		  {
		  int Reject;
		  if (In_Mismatch<6) Reject=TRUE;
		  else//high mismatch situation...
		  {
		  if (Tag.Mismatches > Last_Mismatch_Written)//we already have a better hit..
		  {
		  Hits = 1;
		  return;
		  }
		  if (Gap==1)
		  {
		  Reject=TRUE;//if mismatches are the same as before...
		  if (Tag.Mismatches < Last_Mismatch_Written)//we have a better hit... 
		  {
		  Write_Buffer_Ptr=Last_Write_Buffer;
		  Hits=1;Reject=FALSE;
		  }
		  else if(Tag.Mismatches>Last_Mismatch_Written) {Hits--;return;}//worse hit. Ignore...
		  }
		  else//current hit has many occurances...
		  {
		  Reject=TRUE;
		  }

		  }
		if (Reject)
		{

			Hits=MAXHITS;
			In_Mismatch=7;
			Write_Buffer_Ptr=Last_Write_Buffer;
			return;
		}
	}*/

	New_Record[0]='%';
	New_Record[1]=M.Read.Tag_Number;
	New_Record[2]=Tag.Strand;
	New_Record[3]=M.Read.NCount;
	fwriteX(&New_Record,4,1,Output_File);
	//memcpy(M.Write_Buffer_Ptr,New_Record,4);Write_Buffer_Ptr += 4;
	if (M.Read.NCount)
	{
		fwriteX(&M.Read.N,2*M.Read.NCount,1,Output_File);//;Write_Buffer_Ptr += Skip_Length;//+1;
		//memcpy(Write_Buffer_Ptr,M.Read.N,2*M.Read.NCount);Write_Buffer_Ptr += 2*M.Read.NCount;
	}
	memcpy(MismatchesGIS.Mismatch_Pos,Tag.Mismatch_Pos,Tag.Mismatches);//MAX_MISMATCHES_BOUND);
	MismatchesGIS.Mismatch_Char=Tag.Mismatch_Char;
	Record.Start=Tag.Start;
	Record.Tag=M.Read.Read_Number;
	Record.Skip=Tag.Skip;
	Record.Mismatches=Tag.Mismatches;
	if (ONEFMINDEX) Record.Index=REVERSE; else Record.Index=Tag.FMIndex;
	int Skip_Length=Tag.Mismatches+sizeof(unsigned);
	fwriteX(&Record,sizeof(Record),1,Output_File);
	//memcpy(Write_Buffer_Ptr,&Record,sizeof(Record));Write_Buffer_Ptr += sizeof(Record);
	fwriteX(&MismatchesGIS,Skip_Length,1,Output_File);//;Write_Buffer_Ptr += Skip_Length;//+1;
	//memcpy(Write_Buffer_Ptr,&MismatchesGIS,Skip_Length);Write_Buffer_Ptr += Skip_Length;//+1;
	//Last_Mismatch_Written=Tag.Mismatches;
	}
	funlockfile(Output_File);
}

void Print_Hits(int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char type,Offset_Record *Genome_Offsets,std::string readString,MEMX & M,LEN & L,FILE* Output_File, char ONEFMINDEX,BWT *revfmi,BWT *fwfmi, OUTPUT & O)
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
	
	int i=0;//allHits=M.Hit_Array_Ptr;
	while (i <M.Hit_Array_Ptr)
	{
		SARange Tag=M.Hit_Array[i++];


		if ( ONEFMINDEX && FORWARD==Tag.FMIndex )// forward search index...
		{
			Convert_To_Reverse(Tag,L,revfmi,M);
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
		Hits=Process_GIS(SecondMisN,Multipe,minTrueMis,maxTrueMis,type,readString,O.FILETYPE,fwfmi,revfmi,Record,MismatchesGIS,New_Record,L.STRINGLENGTH,TRUE,O.PLUSSTRAND,O.MaxHits,O.Offset,Output,O.Location_Array,O.Genome_Count,Genome_Offsets,O.Length_Array,M.Read.N,Hits);
	}
	O.Buffer_End=Output;
}


void Convert_To_Reverse(SARange &Tag,LEN & L,BWT *revfmi,MEMX & M)
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
		Search_Forwards_Exact(Current_Tag,Tag,-1,L.STRINGLENGTH,revfmi,M,L.LH);
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
	Location=Location-I->first;
}

int  Location_To_Genome(unsigned & Location,unsigned *Offsets,int Genome_Count)
{
	assert(false);
	int Genome_Position=0;
	while ( Genome_Position< Genome_Count )
	{
		if (Location < Offsets[Genome_Position]) break;
		Genome_Position++;
	}
	Genome_Position--;
	Location=Location-Offsets[Genome_Position];
	return Genome_Position;
}

unsigned Process_GIS(int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char hitType,std::string readString,char FILETYPE,BWT* fwfmi,BWT* revfmi,Output_Record Record,Mismatches_Record_GIS MismatchesGIS,char* New_Record,int StringLength,char USELOCATION,char PLUSSTRAND,unsigned MAXHITS,unsigned Offset,char* & Output,unsigned* Location_Array,int Genome_Count,Offset_Record* Genome_Offsets,int* Length_Array,char* N,unsigned Hits)
{
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	char* Quality[MAXTAG+1];
	unsigned Conversion_Factor;
	char* Mismatch_Desc_Ptr;
	char* Buffer_Index=Output;
	char Mismatches_Desc[1000];
	unsigned Location;
	int Genome_Position;
	Ann_Info A;
	

	int Mismatch_Count = Record.Mismatches;
	StringLength=Length_Array[New_Record[0]];
	Conversion_Factor=revfmi->textLength-StringLength;
	Ins_Format=(char*)"%d<%c\t";Del_Format=(char*)"%d>D\t";Mis_Format=(char*)"%d>%c\t"; 

	//--------------------------- Decode Mismatches ----------------------------------------------
	Mismatch_Desc_Ptr=Mismatches_Desc;

	*Mismatch_Desc_Ptr=0;
	//--------------------------- Decode Mismatches ----------------------------------------------
	
	for (int j=0;j<=Record.Gap && Hits< maxhits-1 ;j++) //MAXHITS //&& Hits< maxhits-1
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
		//Genome_Position=Location_To_Genome(Location,Location_Array,Genome_Count);
		Location_To_Genome(Location,A);
		if(MAX_MISMATCHES>0)
		{
			Mismatch_Count=0;
			int pos =Location+1;
			int H=String_Hash[A.Name];
			char* Genome;
			Genome=Genome_List[H].Genome;//load current genome..
			int Hash_Index=Genome_List[H].Index;//load current genome..
			for(int k=0; k<Length_Array[New_Record[0]]; k++)
			{
				if (pos+k-1 >= Genome_Offsets[Hash_Index+1].Offset) break;
				
				char genome_Char = toupper(Genome[pos+k-1]);//
				
				if (hitType=='1' || hitType=='3') {
					if (!(readString[k]=='C' && genome_Char=='C') && !(readString[k]=='T' && genome_Char=='C') &&  (readString[k] != genome_Char) )
					{
						Mismatch_Count++;
					}
				}
				else if (hitType=='2' || hitType=='4') {
					if (!(readString[k]=='G' && genome_Char=='G') && !(readString[k]=='A' && genome_Char=='G') && (readString[k] != genome_Char) )
					{
						Mismatch_Count++;
					}
				}
			}
			if(Mismatch_Count<minTrueMis) {minTrueMis=Mismatch_Count;Multipe=false;}
			else if(Mismatch_Count==minTrueMis) Multipe=true;
			
			if(Mismatch_Count>maxTrueMis) {maxTrueMis=Mismatch_Count;}
		}
		char CIGr[100];sprintf(CIGr,"%dM",Length_Array[New_Record[0]]);
		if(Mismatch_Count == minTrueMis  ) 
		{
			Output+=sprintf(Output,"%c\t%s\t%c\t%u\t%d\t%d\t%s\n",New_Record[0]+'0',A.Name,New_Record[1],Location+1,Mismatch_Count,Length_Array[New_Record[0]],CIGr);//,Record.Tag,Mismatches_Desc
			Hits++;
		}else if(Mismatch_Count == minTrueMis+1)
		{
			SecondMisN++;
		}
		Record.Start++;
	}

	return Hits;
//--------------------------- Decode Mismatches ----------------------------------------------
}

void Load_Location(char* LOCATIONFILE,std::map <unsigned, Ann_Info> & Annotations)
{
	char Genome_Name_Buf[300];
	char* Genome_Name;unsigned Off_Cum=0,Off=0;int Genome_Count=0;

	FILE* Location_File=File_Open(LOCATIONFILE,"r");
	while (fgets(Genome_Name_Buf,300,Location_File)!=0)
	{
		Off=atoi(Genome_Name_Buf);
		if (Genome_Count)
		{
			Annotations[Off_Cum].Size=Off;
			Annotations[Off_Cum].Name=Genome_Name;
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
	if (!Genome_Count) {fprintf(stderr,"Error Loading genome locations..\n");exit(1);} //else fprintf(stderr,"%d Genomes Loaded..\n",Genome_Count);
}

int Load_Location(char* LOCATIONFILE, Offset_Record* Genome_Offsets,unsigned* Location_Array)
{
	assert(false);
	FILE* Location_File=File_Open(LOCATIONFILE,"r");
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
	Location_Array[0]=0;
	for ( int i=1;i<Genome_Count;i++)
	{
		Location_Array[i]=Location_Array[i-1]+Genome_Offsets[i].Offset;
	}
	return Genome_Count;
}
//}-----------------------------------  Decode ---------------------------------------------------------


//{-----------------------------------  Progress  ---------------------------------------------------------

void Show_Progress(unsigned Percentage)
{
	if (Percentage >=97) return;
	printf("+%u%%\b\b\b",Percentage);
	fflush(stdout);
}
//}-----------------------------------  Progress  ---------------------------------------------------------

