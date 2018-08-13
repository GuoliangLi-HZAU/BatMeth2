//BATMAN 1.1 - added Substring search..
//BATMAN 1.10 - enhanced default output processing...
//		handle N's.
//		print blanks...
//		--maxhits bug fixes
//{-----------------------------  INCLUDE FILES  -------------------------------------------------/
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include "bfix.h"
//#include <dvec.h>
#include <getopt.h>
extern "C" 
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}
//}-----------------------------  INCLUDE FILES  -------------------------------------------------/

//{-----------------------------  DEFINES  -------------------------------------------------/
//#define BUILTIN_LOG
#define DEBUG
#define SAGAP_CUTOFF 2
#define FIELD_LENGTH 32
#define MAXGAP 30000 //maximum SA range to sort...
#define DEFAULT 0
#define DEEP 1
#define INTEGERSIZE 8 //integer size
#define FALSE 0
#define TRUE 1
//}-----------------------------  DEFINES  -------------------------------------------------/

//{-----------------------------  STRUCTS  -------------------------------------------------/
struct SARange
{

	unsigned Start;
	unsigned End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
	int Tag;//Tag number
	unsigned Mismatch_Pos;//6|6|...
	unsigned Mismatch_Char;//2|2|...
	
};

struct Range
{
	unsigned Start;
	unsigned End;
	int Label;//Final Label of the range...
};

struct Mismatches_Record
{

	int 	Gap;
	int Mismatch_Pos;//6|6|...
	unsigned Mismatch_Char;//2|2|...


}__attribute__((__packed__));

#define REVERSE 1
#define FORWARD 0


struct Branches
{
	char  Is_Branch [4];
};

struct SA
{
	unsigned Start;
	unsigned End;
	unsigned Start_Location;// exact location of the first occurance...
	unsigned End_Location;
} Suffix_Range;

struct Range_Record
{
	int Offset;
	unsigned Location;
}; 

//}-----------------------------  STRUCTS  -------------------------------------------------/*

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

void Build_Tables();
void Load_Indexes();
void Open_Files();
void Build_Indices();
void Print_Indices();
void Build_Preindex_Forward(Range Range, int Level, int Bit);
void Build_Preindex_Backward(Range Range, int Level, int Bit);
void Parse_Command_line(int argc, char* argv[]);
void Show_Progress(float Percentage);
int SA_cmp(const void *First,const void *Second);
void fwriteX(void* Offset,unsigned long Size1,unsigned long Size2,FILE* Handle);
int Range_Record_cmp(const void *First,const void *Second);
unsigned log2(unsigned v); // 32-bit word to find the log of
FILE* File_Open(const char* File_Name,const char* Mode);
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName, const char* SAFile); 
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

//{---------------------------- GLOBAL VARIABLES -------------------------------------------------

FILE* SaRanges;
FILE* Ranges;
FILE* Index;
FILE* Blocks;
FILE* Log_SFile;

char* LOG_SUCCESS_FILE=NULL;

MMPool *mmPool;
BWT *fwfmi,*revfmi;
off64_t File_Size;
time_t Start_Time,End_Time;

SA* SA_Array;
Range_Record Range_Array[MAXGAP];
unsigned Sorted_SA_Array[MAXGAP];

int Actual_Tag;
int LOOKUPSIZE=6;
int HITMODE = DEFAULT;

unsigned Hits,Total_Hits=0;
unsigned Forward_Start_Lookup[4],Forward_End_Lookup[4],Backward_Start_Lookup[4],Backward_End_Lookup[4];
unsigned* Forward_Start_LookupX;
unsigned* Forward_End_LookupX;
unsigned* Backward_Start_LookupX;
unsigned* Backward_End_LookupX;
unsigned* Print_Stack;
unsigned DISKBUFFERSIZE =1000;

char BWTFILE_DEFAULT[] = "genome.bwt"; 
char OCCFILE_DEFAULT[] ="genome.fmv";
char REVBWTINDEX_DEFAULT[] ="revgenome.bwt";
char REVOCCFILE_DEFAULT[] ="revgenome.fmv";
char GENOMEFILE_DEFAULT[]="genome";
char OUTFILE_DEFAULT[]="hits.txt";
char SKIP_SA_ENUM=FALSE;
char COMPRESS=FALSE;

char* Source;
char* BWTFILE ; 
char* OCCFILE ;
char* SAFILE ;
char* REVSAFILE ;
char* REVBWTINDEX;
char* REVOCCFILE;
char* GENOMEFILE;
char* OUTFILE;
char* BLKFILE;
char* INDFILE;
char* RANGEFILE;
char* SORTEDRANGEFILE;
 
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------

//{---------------------------- Command Line  -------------------------------------------------
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"output",1,NULL,'o'},
{"genome",1,NULL,'g'},
{"stringlength",1,NULL,'l'},
{"buffersize",1,NULL,'b'},
{"skipenum",0,NULL,'e'},
{"compress",0,NULL,'c'},
{0,0,0,0}
};

//}---------------------------- Command Line -------------------------------------------------
int main(int argc, char* argv[])
{
	
	time(&Start_Time);
	LOOKUPSIZE=18;
	Parse_Command_line(argc,argv);	
	Load_Indexes();	
	printf("Unique hit length %d ..\n",LOOKUPSIZE);
	//Hits=0;Total_Hits=0;
	//printf("===============================================================]\r[");	
	//fflush(stdout);

	//Ranges=File_Open("RangesF","wb");
	Ranges=File_Open(SORTEDRANGEFILE,"wb");
	//Index=File_Open("indexF.dat","wb");
	Index=File_Open(INDFILE,"wb");
	//Blocks=File_Open("blocksF.dat","wb");
	Blocks=File_Open(BLKFILE,"wb");
	if(!SKIP_SA_ENUM)
	{
		SaRanges=File_Open(RANGEFILE,"wb");
		printf("===============================================================]\r[");	
		// Enumerate SA Ranges..
		Build_Tables();
		fwriteX(&Hits,1,sizeof(Hits),SaRanges);
		fclose(SaRanges);
		printf("\r[+++++++++++++++++++++++++++ 100%% ++++++++++++++++++++++++++++++]\n");	
		time(&End_Time);printf("\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}

//Build the final indices...
	time(&Start_Time);
	Build_Indices();
	//Print_Indices();
	printf("\r[++++++++100%%+++++++++]\n");//progress bar....
	time(&End_Time);printf("\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	
	
}

/*void Build_Indices()
{

	SaRanges=File_Open("SA-F","rb");
	fseek(SaRanges,-sizeof(Hits),SEEK_END);
	fread(&Hits,sizeof(Hits),1,SaRanges);
	SA_Array=(SA*)malloc(sizeof(SA)*Hits);
	fseek(SaRanges,0,SEEK_SET);//go top
	fread(&SA_Array[0],sizeof(SA),Hits,SaRanges);
	qsort(&SA_Array[0], Hits, sizeof(SA), SA_cmp);
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Ranges);//write Ranges structure...

	unsigned long Block_File_Size=0;
	unsigned Progress;
	unsigned Percent_Mark=Hits/20;
	int Gap;
	printf("%d SA Ranges indexed ...\n",Hits);
	printf("======================]\r[");//progress bar....

	for (unsigned i=0;i<Hits;i++)
	{
		Gap = SA_Array[i].End-SA_Array[i].Start+1;
		//if (Gap<MAXGAP) Block_File_Size += log2(Gap)*(Gap-1);
		Block_File_Size += Gap;
	}
	Block_File_Size=Block_File_Size/8+1;
	unsigned Conversion_Factor=revfmi->textLength-LOOKUPSIZE;

	unsigned char* Block_Structure=(unsigned char*)malloc(Block_File_Size);//66550747); 
	if (NULL==Block_Structure) {printf("malloc error:\n");exit(1);}
	unsigned Block_Index=1;
	int Field_Length=32;
	for (unsigned i=0;i<Hits;i++)
	{
		Progress++;
		if(Progress > Percent_Mark) {Show_Progress(i*100/Hits);Progress=0;};
		int Gap = SA_Array[i].End-SA_Array[i].Start+1;
		if (Gap <MAXGAP )
		{ 
			//	printf("_______________________\n");
			for (int j=0;j<Gap;j++) 
			{
				Range_Array[j].Offset=j;
				Range_Array[j].Location=Conversion_Factor-BWTSaValue(revfmi,SA_Array[i].Start+j);
				//printf("%d\n", Range_Array[j].Location);
			}
			qsort(&Range_Array[0], Gap, sizeof(Range_Record), Range_Record_cmp);
			SA_Array[i].End=Block_Index;
			SA_Array[i].Start_Location=Range_Array[0].Location;
			//printf("%d:",SA_Array[i].Start_Location);
			//bfi(Block_Structure,Block_Index,FIELD_LENGTH,Range_Array[0].Location);Block_Index=Block_Index+FIELD_LENGTH;
			Field_Length=log2(Gap);
			//Block_Structure[Block_Index]=Range_Array[0].Location;
			//Block_Index++;
			for (int j=1;j<Gap-1;j++) 
			{
				//Block_Structure[Block_Index]=Range_Array[j].Location;
				//Block_Index++;
				bfi(Block_Structure,Block_Index,Field_Length,Range_Array[j].Offset);
				Block_Index=Block_Index+Field_Length;
				//printf("%d:",Range_Array[j].Location);
			}
			SA_Array[i].End_Location=Range_Array[Gap-1].Location;
			//printf("%d:",SA_Array[i].End_Location);
			//	printf("_______________________\n");
			//printf("%d : %d\n", SA_Array[i].Start,SA_Array[i].End);
		}
	}

	//printf("\r[++++++++100%%+++++++++]\n");//progress bar....
	fwriteX(&Hits,1,sizeof(unsigned),Index);
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Index);//write Ranges structure...
	fwriteX(Block_Structure,1,Block_File_Size,Blocks);
}*/

void Build_Indices()
{

	SaRanges=File_Open(RANGEFILE,"rb");
	fseek(SaRanges,-sizeof(Hits),SEEK_END);
	fread(&Hits,sizeof(Hits),1,SaRanges);
	SA_Array=(SA*)malloc(sizeof(SA)*Hits);
	fseek(SaRanges,0,SEEK_SET);//go top
	fread(&SA_Array[0],sizeof(SA),Hits,SaRanges);
	qsort(&SA_Array[0], Hits, sizeof(SA), SA_cmp);
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Ranges);//write Ranges structure...
	fwriteX(&COMPRESS,1,1,Blocks);//write header...
	fwriteX(&Hits,1,sizeof(Hits),Blocks);//write header...

	unsigned long Block_File_Size=0;
	unsigned Progress;
	unsigned Percent_Mark=Hits/20;
	int Gap;
	if (COMPRESS) printf("Building Compressed index...\n"); else printf ("Building uncompressed index...\n");
	printf("%d SA Ranges indexed ...\n",Hits);
	printf("======================]\r[");//progress bar....

	for (unsigned i=0;i<Hits;i++)
	{
		Gap = SA_Array[i].End-SA_Array[i].Start+1;
		if (COMPRESS)
		{
			if (Gap<MAXGAP) Block_File_Size += log2(Gap)*(Gap-1);
		}
		else
			Block_File_Size += Gap;
	}
	if (COMPRESS) Block_File_Size=Block_File_Size/8+1;
	unsigned Conversion_Factor=revfmi->textLength-LOOKUPSIZE;

	unsigned char* Block_Structure=(unsigned char*)malloc(Block_File_Size);//66550747); 
	if (NULL==Block_Structure) {printf("malloc error:\n");exit(1);}
	unsigned Block_Index;
	if(COMPRESS) Block_Index=1; else Block_Index=0;
	int Field_Length=32;
	for (unsigned i=0;i<Hits;i++)
	{
		Progress++;
		if(Progress > Percent_Mark) {Show_Progress(i*100/Hits);Progress=0;};
		int Gap = SA_Array[i].End-SA_Array[i].Start+1;
		if (Gap <MAXGAP )
		{ 
			for (int j=0;j<Gap;j++) 
			{
				Range_Array[j].Offset=j;
				Range_Array[j].Location=Conversion_Factor-BWTSaValue(revfmi,SA_Array[i].Start+j);
			}
			qsort(&Range_Array[0], Gap, sizeof(Range_Record), Range_Record_cmp);
			SA_Array[i].End=Block_Index;
			SA_Array[i].Start_Location=Range_Array[0].Location;
			Field_Length=log2(Gap); if (Field_Length==32) printf("Error");
			int k=0;
			for (int j=1;j<Gap-1;j++) 
			{
				if (COMPRESS)
				{
					bfi(Block_Structure,Block_Index,Field_Length,Range_Array[j].Offset);
					Block_Index=Block_Index+Field_Length;
				}
				else
				{
					Sorted_SA_Array[j-1]=Range_Array[j].Location;
					Block_Index++;k++;
				}
			}
			if (!COMPRESS) fwriteX(Sorted_SA_Array,k,sizeof(unsigned),Blocks);

			SA_Array[i].End_Location=Range_Array[Gap-1].Location;
		}
	}

	//printf("\r[++++++++100%%+++++++++]\n");//progress bar....
	////fwriteX(&Hits,1,sizeof(unsigned),Index);
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Index);//write Ranges structure...
	if(COMPRESS) fwriteX(Block_Structure,1,Block_File_Size,Blocks);
}

/*void Print_Indices()
{

	SaRanges=File_Open("SA-F","rb");
	FILE* Prindex=File_Open("prindex.txt","w");
	fseek(SaRanges,-sizeof(Hits),SEEK_END);
	fread(&Hits,sizeof(Hits),1,SaRanges);
	SA_Array=(SA*)malloc(sizeof(SA)*Hits);
	fseek(SaRanges,0,SEEK_SET);//go top
	fread(&SA_Array[0],sizeof(SA),Hits,SaRanges);
	qsort(&SA_Array[0], Hits, sizeof(SA), SA_cmp);
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Ranges);//write Ranges structure...

	unsigned long Block_File_Size=0;
	unsigned Progress;
	unsigned Percent_Mark=Hits/20;
	int Gap;
	printf("%d SA Ranges indexed ...\n",Hits);
	printf("======================]\r[");//progress bar....

	for (int i=0;i<Hits;i++)
	{
		Gap = SA_Array[i].End-SA_Array[i].Start+1;
		if (Gap<MAXGAP) Block_File_Size += log2(Gap)*(Gap-1);
	}
	Block_File_Size=Block_File_Size/8+1;
	unsigned Conversion_Factor=revfmi->textLength-LOOKUPSIZE;

	unsigned char* Block_Structure=(unsigned char*)malloc(Block_File_Size);//66550747); 
	if (NULL==Block_Structure) {printf("malloc error:\n");exit(1);}
	unsigned Block_Index=1;
	int Field_Length=32;
	for (int i=0;i<Hits;i++)
	{
		Progress++;
		if(Progress > Percent_Mark) {Show_Progress(i*100/Hits);Progress=0;};
		int Gap = SA_Array[i].End-SA_Array[i].Start+1;
		if (Gap <MAXGAP )
		{ 
			fprintf(Prindex,"_______________________\n");
			for (int j=0;j<Gap;j++) 
			{
				Range_Array[j].Offset=j;
				Range_Array[j].Location=Conversion_Factor-BWTSaValue(revfmi,SA_Array[i].Start+j);
				fprintf(Prindex,"%d\n", Range_Array[j].Location);
			}
			qsort(&Range_Array[0], Gap, sizeof(Range_Record), Range_Record_cmp);
			SA_Array[i].End=Block_Index;
			SA_Array[i].Start_Location=Range_Array[0].Location;
			fprintf(Prindex,"%d:",SA_Array[i].Start_Location);
			//bfi(Block_Structure,Block_Index,FIELD_LENGTH,Range_Array[0].Location);Block_Index=Block_Index+FIELD_LENGTH;
			Field_Length=log2(Gap);
			//Block_Structure[Block_Index]=Range_Array[0].Location;
			//Block_Index++;
			for (int j=1;j<Gap-1;j++) 
			{
				//Block_Structure[Block_Index]=Range_Array[j].Location;
				//Block_Index++;
				bfi(Block_Structure,Block_Index,Field_Length,Range_Array[j].Offset);
				Block_Index=Block_Index+Field_Length;
				fprintf(Prindex,"%d:",Range_Array[j].Location);
			}
			SA_Array[i].End_Location=Range_Array[Gap-1].Location;
			fprintf(Prindex,"%d:",SA_Array[i].End_Location);
			fprintf(Prindex,"_______________________\n");
			fprintf(Prindex,"%d : %d\n", SA_Array[i].Start,SA_Array[i].End);
		}
	}

	//printf("\r[++++++++100%%+++++++++]\n");//progress bar....
	fwriteX(&SA_Array[0],Hits,sizeof(SA),Index);//write Ranges structure...
	fwriteX(Block_Structure,1,Block_File_Size,Blocks);
}*/
static const char LogTable256[] = 
{

  1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  
};
 
inline unsigned log2(unsigned v) // 32-bit word to find the log of
{
#ifdef BUILTIN_LOG
	return (32 -__builtin_clz(v));
#else
	unsigned r;     // r will be lg(v)
	register unsigned int t, tt; // temporaries

	if (tt = v >> 16)
	{
		r= (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	}
	else 
	{
		r= (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
	}
	return r;
#endif
}

int SA_cmp(const void *First,const void *Second)
{
	SA* F=(SA*) First;
	SA* L=(SA*) Second;
	if (F->Start > L->Start) return 1; else return -1;
	//return (*F).Start-L->Start;
	/* integer comparison: returns negative if b > a 
	and positive if a > b */
}

int Range_Record_cmp(const void *First,const void *Second)
{
	Range_Record* F=(Range_Record*) First;
	Range_Record* L=(Range_Record*) Second;
	if (F->Location > L->Location) return 1; else return -1;
	//return F->Location-L->Location;
	/* integer comparison: returns negative if b > a 
	and positive if a > b */
}


void Build_Preindex_Backward(Range Range, int Level, int Bit)
{
	
	if (LOOKUPSIZE==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		//Backward_Start_LookupX[Range.Label] = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
		//Backward_End_LookupX[Range.Label] = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
		unsigned S=fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
		unsigned L=fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
		if(L>S)//(L>=S)  
		{
			Total_Hits++;
			if(L-S<=SAGAP_CUTOFF)
			{
			}
			else
			{
				Suffix_Range.Start=S;Suffix_Range.End=L;
				fwriteX(&Suffix_Range,1,sizeof(Suffix_Range),SaRanges);
				Hits++;
			}
		}
		//if(S<=L) printf("%d : %d\n",S,L);
	}
	else
	{


		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		Range.Start = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.Start , Bit) + 1;
		Range.End = fwfmi->cumulativeFreq[Bit] + BWTOccValue(fwfmi, Range.End+1, Bit);
		if (Range.End > Range.Start)//(Range.End >= Range.Start)
		{
			Level ++;
			for ( int i=0;i<4;i++)
			{
				Build_Preindex_Backward( Range, Level,i);
			}
		}
	}
}

void Build_Preindex_Forward(Range Range, int Level, int Bit)
{
	static float Percentage=0;
	int debug;
	if (LOOKUPSIZE==Level) 
	{
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label
		unsigned S= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
		unsigned L= revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
		if(L>S)//if(L>=S) We keep only non-unique hits.. 
		{
			Total_Hits++;
			if(L-S<=SAGAP_CUTOFF)
			{
				debug=1;
			}
			else
			{
				Suffix_Range.Start=S;Suffix_Range.End=L;
				fwriteX(&Suffix_Range,1,sizeof(Suffix_Range),SaRanges);
				Hits++;
			}
		}
	}
	else
	{

		if (Level==3) {Percentage += 1.5;Show_Progress(Percentage);}
		Range.Label=Range.Label | (Bit<<2*(Level-1));//Calculate label 
		Range.Start = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.Start , Bit) + 1;
		Range.End = revfmi->cumulativeFreq[Bit] + BWTOccValue(revfmi, Range.End+1, Bit);
		if (Range.End > Range.Start)//(Range.End >= Range.Start) only non unique needed...
		{
			Level ++;
			for ( int i=0;i<4;i++)
			{
				Build_Preindex_Forward( Range, Level,i);
			}
		}
	}

}
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
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

//}----------------------------------- FILE HANDLING ---------------------------------------------------------

//{----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initFMI
 *  Description:  Opens FM index fmiFile
 * =====================================================================================
 */
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName, const char* SAFile) 

{
	BWT *fmi;
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	fmi = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SAFile, NULL, NULL, NULL);//Load FM index

	return fmi;
}

//}----------------------------------- FM INDEX ROUTINES ---------------------------------------------------------

//{-----------------------------  Parse Command Line  -------------------------------------------------
void Parse_Command_line(int argc, char* argv[])
{
	int Current_Option=0;
	char* Short_Options ="ho:b:g:l:ec";//allowed options....
	char* This_Program = argv[0];//Current program name....
	char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --output | -o <filename>\t\t Name of output file\n"
" --genome | -g <filename>\t\t Name of the reference genome\n"
" --stringlength | -l <integer> \t\t Length of strings to build index for..\n"
" --buffersize | -b <integer> \t\t Size of disk buffers\n"
" --skipenum | -e \t\t\t Skip enumeration of SA ranges..\n"
" --compress | -c \t\t\t Compress index..\n"
;

	Source=(char*)malloc(sizeof(char)*6000);//create space for file names...
	char *options, *value; 
	char* Name;int Last_Dash;char *Ind,*Blk;char* Genome_Name;

	OUTFILE=OUTFILE_DEFAULT;GENOMEFILE=GENOMEFILE_DEFAULT;
	BWTFILE = BWTFILE_DEFAULT; 
	OCCFILE = OCCFILE_DEFAULT;
	REVBWTINDEX=REVBWTINDEX_DEFAULT;
	REVOCCFILE=REVOCCFILE_DEFAULT;

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'l':
				LOOKUPSIZE=atoi(optarg);
				break;
			case 'e':
				SKIP_SA_ENUM = TRUE;
				break;
			case 'o':
				OUTFILE=optarg;
				break;
			case 'b':
				DISKBUFFERSIZE=atol(optarg);
				break;
			case 'c':
				COMPRESS=TRUE;
				break;
			case 'g':
				Name=optarg;Last_Dash=0;Genome_Name=optarg;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

				REVBWTINDEX = (char*)Source;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(REVBWTINDEX,optarg,Last_Dash);
				REVBWTINDEX[Last_Dash+0]='r';REVBWTINDEX[Last_Dash+1]='e';REVBWTINDEX[Last_Dash+2]='v';
				strcpy(REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(REVBWTINDEX+Last_Dash+3,".bwt"); 

				BWTFILE=REVBWTINDEX+600;
				strncpy(BWTFILE,optarg,Last_Dash);
				strcpy(BWTFILE+Last_Dash,Genome_Name+1);
				strcat(BWTFILE+Last_Dash,".bwt"); 


				REVOCCFILE = BWTFILE+600;
				strncpy(REVOCCFILE,optarg,Last_Dash);
				REVOCCFILE[Last_Dash+0]='r';REVOCCFILE[Last_Dash+1]='e';REVOCCFILE[Last_Dash+2]='v';
				strcpy(REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVOCCFILE+Last_Dash+3,".fmv"); 


				OCCFILE=REVOCCFILE+600;			
				strncpy(OCCFILE,optarg,Last_Dash);
				strcpy(OCCFILE+Last_Dash,Genome_Name+1);
				strcat(OCCFILE+Last_Dash,".fmv"); 

				SAFILE=OCCFILE+500;			
				strncpy(SAFILE,optarg,Last_Dash);
				strcpy(SAFILE+Last_Dash,Genome_Name+1);
				strcat(SAFILE+Last_Dash,".sa");

				REVSAFILE = SAFILE+500;
				strncpy(REVSAFILE,optarg,Last_Dash);
				REVSAFILE[Last_Dash+0]='r';REVSAFILE[Last_Dash+1]='e';REVSAFILE[Last_Dash+2]='v';
				strcpy(REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(REVSAFILE+Last_Dash+3,".sa"); 

				BLKFILE = REVSAFILE+500;
				strncpy(BLKFILE,optarg,Last_Dash);
				strcpy(BLKFILE+Last_Dash,Genome_Name+1);
				strcat(BLKFILE+Last_Dash,".blk."); 
				Blk=BLKFILE+strlen(BLKFILE);

				INDFILE = BLKFILE+500;
				strncpy(INDFILE,optarg,Last_Dash);
				strcpy(INDFILE+Last_Dash,Genome_Name+1);
				strcat(INDFILE+Last_Dash,".ind."); 
				Ind=INDFILE+strlen(INDFILE);

				RANGEFILE = INDFILE+500;
				strncpy(RANGEFILE,optarg,Last_Dash);
				strcpy(RANGEFILE+Last_Dash,Genome_Name+1);
				strcat(RANGEFILE+Last_Dash,".range"); 

				SORTEDRANGEFILE = RANGEFILE+500;
				strncpy(SORTEDRANGEFILE,optarg,Last_Dash);
				strcpy(SORTEDRANGEFILE+Last_Dash,Genome_Name+1);
				strcat(SORTEDRANGEFILE+Last_Dash,".sort"); 

				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
	if (LOOKUPSIZE>0)
	{
		sprintf(Ind,"%d",LOOKUPSIZE);
		sprintf(Blk,"%d",LOOKUPSIZE);
	}
	else
	{
		printf("Parse_Command_line():Bad String Length\n");
	}
}

//}-----------------------------  Parse Command Line  -------------------------------------------------
void Build_Tables()
{
	Range LRange;
	for( int i=0;i<4;i++)//fill lookup tables for first character...
	{
		/*Forward_Start_Lookup[i]=revfmi->cumulativeFreq[i] + 1;
		Forward_End_Lookup[i]=revfmi->cumulativeFreq[i + 1];
		Backward_Start_Lookup[i]=fwfmi->cumulativeFreq[i] + 1;
		Backward_End_Lookup[i]=fwfmi->cumulativeFreq[i + 1];//
*/
		LRange.Start=revfmi->cumulativeFreq[i] + 1;
		LRange.End=revfmi->cumulativeFreq[i + 1];
		LRange.Label=i;
		if (LRange.Start <=LRange.End) for(int j=0;j<4;j++) Build_Preindex_Forward(LRange, 2, j);

		/*LRange.Start=fwfmi->cumulativeFreq[i] + 1;
		LRange.End=fwfmi->cumulativeFreq[i + 1];
		LRange.Label=i;
		for(int j=0;j<4;j++) Build_Preindex_Backward(LRange, 2, j);*/

	}
}

void Load_Indexes()
{
	fwfmi=initFMI(BWTFILE,OCCFILE,SAFILE);//Load FM indexes
	revfmi=initFMI(REVBWTINDEX,REVOCCFILE,REVSAFILE);
	unsigned SOURCELENGTH = fwfmi->textLength;
	if (SOURCELENGTH!=revfmi->textLength)
	{ 
		printf("FM index load error \n"); 
		exit(1);
	}
	//FWDInverseSA0=fwfmi->inverseSa0;
}

void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}


void fwriteX(void* Offset,unsigned long Size1,unsigned long Size2,FILE* Handle)
{
		if(Size2 !=fwrite(Offset,Size1,Size2,Handle)) {printf("Error writing to file...\n");exit(1);};
}
