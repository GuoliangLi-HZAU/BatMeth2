#ifndef BATFUNC
#define BATFUNC

#include <limits.h>
#include <string>
#include <string.h>
#include <map>
#include "zlib.h"
#include "common.h"
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}

#define QUIT_TAG INT_MAX
struct Ann_Info
{
	unsigned Size;
	unsigned Cumulative_Size;
	char* Name;
};
void Initial_Length(LEN & L);
void Analyze_File(INFILE & I);
void GET_LEN(LEN & L,READ & Read);
void Analyze_File(INFILE & I,LEN & L);
void Write_Header(HEADER & Header,FILE* Output_File,BATPARAMETERS P,INFILE F,LEN L);
size_t fwriteX ( const void * ptr, size_t size, size_t count, void* stream );
void fprintfX(void* Handle,char* Format, char* String);

int Map_Strand(int Mismatches,int MAXHITS,LEN & L,BWT *fwfmi,BWT* revfmi,MEMX & M);
unsigned Map_Strand_Guess(int MAX_MISMATCHES,int MAXHITS,LEN & L,BWT *fwfmi,BWT* revfmi,MEMX & M);

void Load_Indexes(BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,FMFILES & FM);//char *BWT,char *OCC,char *SA,char *REVBWT,char *REVOCC,char *REVSA)
#define Load_FM(fwfmix,revfmix,mmPoolx,FMFilesx) do {FMFilesx.SAFILE=FMFilesx.REVSAFILE=NULL;Load_Indexes(fwfmix,revfmix,mmPoolx,FMFilesx);} while(0)
#define Load_FM_and_Sample(fwfmix,revfmix,mmPoolx,FMFilesx) Load_Indexes(fwfmix,revfmix,mmPoolx,FMFilesx)
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName,const char* SAFile, MMPool *mmPool);
SARange Get_SARange( char New_Char,struct SARange Range,BWT *fmi);
void Get_SARange_Fast( char New_Char,struct SARange & Range,BWT *fmi);
void Get_SARange_Fast_2( char New_Char,struct SARange & Start_Range, struct SARange & Dest_Range,BWT *fmi);
void Split_Read(const int STRINGLENGTH, LEN & L); 
void Init(char Solid,char Npolicy);
void Allocate_Memory(MEMX & M);
void Copy_MEM(MEMLOOK & M,MEMX & F,MEMX & C,char MAX_MISMATCHES);
void Build_Tables(BWT *fwfmi, BWT *revfmi,MEMLOOK & M);
void Build_Preindex_Backward(BWT *fwfmi,Range Range, int Level, int Bit, const MEMLOOK & M);
void Build_Preindex_Forward(BWT *revfmi, Range Range, int Level, int Bit,  const MEMLOOK & M);

char Read_Tag(FILE *Input_File,const char FILETYPE,READ & Read );
bool Process_Read(READ & R, BATREAD & B,MEMX & MF, MEMX & MC);

FILE* File_Open(const char* File_Name,const char* Mode);
void File_OpenZ(const char* File_Name,const char* Mode,gzFile & Handle);
size_t fwriteX ( const void * ptr, size_t size, size_t count, void* stream );
void fprintfX(void* Handle,char* Format, char* String);

unsigned Guess_Orientation(BWT* fwfmi,BWT* revfmi,MEMX & MF,MEMX & MC,LEN & L,GUESS & G,BATREAD & B);
unsigned Zero_Mismatch(char* Current_Tag,LEN L, BWT *revfmi,MEMX & M);
unsigned One_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M);
unsigned Two_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M);
unsigned Three_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M);
unsigned Four_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M);
unsigned Five_Mismatch(char* Current_Tag,LEN L, int MAXHITS, BWT* fwfmi, BWT* revfmi,MEMX & M);

void Search_Backwards(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,MEMX & M);
void Search_Backwards_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M);
void Search_Forwards(const char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT* & revfmi,MEMX & M);
void Search_Exact(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Forwards_OneSA(const char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,MEMX & M);

void Search_Forwards_Exact(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi,MEMX & M,int LH);
void Search_Backwards_Exact(const char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi,MEMX & M);
void Search_Backwards_Exact_X0(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Forwards_0X(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi);
void Search_Backwards_X10(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fmi,MEMX & M);
void Search_Backwards_X10_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fmi,MEMX & M);
void Search_X01(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT* fwfmi,BWT *revfmi,MEMX & M);
void Search_X01_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,LEN & L,BWT *fwfmi,BWT *revfmi,MEMX & M);
void Search_01X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M);
void Search_01X_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M);
void Search_10X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L, int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M);
void Search_10X_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,BWT *revfmi,MEMX & M);
void Search_11X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *revfmi,MEMX & M);
void Search_Half_Tag_11X(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M);
void Search_X11(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,unsigned MAXHITS, LEN & L,BWT *fmi,MEMX & M);
void Search_Half_Tag_X11(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M);
void Search_Half_Tag_X11_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,LEN & L,int MAXHITS,BWT *fmi,MEMX & M);
void Search_01LX(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *revfmi,LEN & L,MEMX & M);
void Search_01LX_OneSA(char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fmi,LEN & L,MEMX & M);
void Search_10LX(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M);
void Search_Backwards_XL10(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *fmi,LEN & L,MEMX & M);
void Search_Backwards_XL10(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS, BWT *fmi,LEN & L,MEMX & M);
void Search_XL01(char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,int MAXHITS,BWT *fwfmi,BWT *revfmi,LEN & L,MEMX & M);

void Reverse(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,LEN & L, BWT *fwfmi,MEMX & M);
void Backwards(char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *revfmi,MEMX & M);
void Branch_Detect (const struct SARange Tag,BWT *fmi,int Start,unsigned* Branch_Characters,SARange* Branch_Ranges);
void Branch_Detect_Backwards (const char* Current_Tag,const struct SARange Tag,BWT *fmi,int Start,unsigned* Branch_Characters,SARange *Branch_Ranges);
void Print_LocationX(SARange & Tag, MEMX & M);
void Convert_To_Reverse(SARange &Tag,LEN & L,BWT *revfmi,MEMX & M);
void Print_Hits(int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char type,Offset_Record *Genome_Offsets,std::string readString,MEMX & M,LEN & L,FILE* Output_File, char ONEFMINDEX,BWT *revfmi,BWT *fwfmi, OUTPUT & O);
void Print_Hits_To_BAT(MEMX & M,LEN & L,FILE* Output_File, char ONEFMINDEX,BWT *revfmi,char FILETYPE);

int Get_Lookup_Size(char MAX_MISMATCHES,char STRINGLENGTH);

//{---------------------------------------- DECODE ROUTINES ---------------------------------------------
void Open_BAT_File(char* INPUTFILE, FILE* & Data_File,In_File & I);
int Load_Location(char* LOCATIONFILE, Offset_Record* Genome_Offsets,unsigned* Location_Array);
void Load_Location(char* LOCATIONFILE,std::map <unsigned, Ann_Info> & Annotations);
int  Location_To_Genome(unsigned & Location,unsigned *Offsets,int Genome_Count);
void Location_To_Genome(unsigned & Location,Ann_Info & A);
unsigned Process_GIS(int & SecondMisN,bool &Multipe,int &minTrueMis,int &maxTrueMis,char hitType,std::string readString,char FILETYPE,BWT* fwfmi,BWT* revfmi,Output_Record Record,Mismatches_Record_GIS MismatchesGIS,char* New_Record,int StringLength,char USELOCATION,char PLUSSTRAND,unsigned MAXHITS,unsigned Offset,char* & Output,unsigned* Location_Array,int Genome_Count,Offset_Record* Genome_Offsets,int* Length_Array,char* N,unsigned Hits);
//}---------------------------------------- DECODE ROUTINES ---------------------------------------------

//{---------------------------------------- MISC ---------------------------------------------
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....
void Show_Progress(unsigned Percentage);
//}---------------------------------------- MISC ---------------------------------------------

#endif

