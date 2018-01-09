#ifndef UTILS
#include "common.h"
//{-----------------------------  DEFINES  -------------------------------------------------/
# define MAX_HITS_TO_STORE 200
#define MAX_MISMATCHES_BOUND 16 //upper boud for mismatch number....
#define TAB	1
#define FQ	2
#define FA	3
#define TWOFILE	4
#define GISMODE 100
#define BUILTIN_LOG
#define UTILS
#define DEBUG
#define SAGAP_CUTOFF 2
#define MAXCOUNT_DEFAULT 1 //30000
#define FIELD_LENGTH 32
#define MAXGAP 30000 //maximum SA range to sort...
#define DEFAULT 0
#define DEEP 1
#define PAIREND 2
#define PAIR_END_SEPERATOR '\t'
#define MAXDES 500 
//#define MAXTAG 180
#define EXTRA 0//Stuff like cr/lf at string seperators in the input file
#define INTEGERSIZE 8 //integer size
#define PACKEDBITSIZE 2
#define BITMASK 1 //3
#define FALSE 0
#define TRUE 1
#define NOMISMATCHES 100 
#define PRINTSTACKSIZE 50
#define MAXSTRINGLENGTH 180 //36//36//6+EXTRA//6+EXTRA//36
#define BUILTIN_LOG

#define START_OF_MARK 9//18  //Start of scanning the tag for the pruning of one mismatch...
//#define BRANCHTHRESHOLD 80 //30 //Threshold at which to check the BWT instead of branching
const int SW_MAX_MIS_TO_STORE=10;
//}-----------------------------  DEFINES  -------------------------------------------------/

//{-----------------------------  STRUCTS  -------------------------------------------------/



struct Record_Info
{
	char Description[MAXDES+1];
	char Quality[MAXDES+1];
	char Positive_Head[MAXDES+1];
	char Reverse_Head[MAXDES+1];
	char Positive_Tail[MAXDES+1];
	char Reverse_Tail[MAXDES+1];
};

struct Tag_Info
{
	unsigned SA_Start;//start of the Sa range
	unsigned Gap;//length of SA range
	unsigned Block_Start;//Start of block info..
	unsigned Index;//Index to SA_index
	unsigned First;//first location of hit
	unsigned Last;//last location of hit
	unsigned Field_Length;
};

struct Pair
{
	unsigned Head;
	unsigned Tail;
	char MismatchesH;
	char MismatchesT;
	unsigned char Bad_Loc[SW_MAX_MIS_TO_STORE];
	char Bad_Char[10];
	char Edit_Distance;
};

#define REVERSE 1
#define FORWARD 0



struct SA
{

	unsigned Start;
	unsigned End;
	unsigned Start_Location;// exact location of the first occurance...
	unsigned End_Location;

};

struct Range_Record
{
	int Offset;
	unsigned Location;
}; 

//}-----------------------------  STRUCTS  -------------------------------------------------/*
#include <stdlib.h>
extern "C"
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}



unsigned Get_File_Size(FILE* File);
unsigned log2(unsigned v); // 32-bit word to find the log of
FILE* File_Open(const char* File_Name,const char* Mode);
BWT* initFMI(const char* BWTCodeFileName,const char* BWTOccValueFileName, const char* SAFile);
SARange Get_SARange( long New_Char,struct SARange Range,BWT *fmi);
void Get_SARange_Fast( long New_Char,struct SARange & Range,BWT *fmi);
void Get_SARange_Fast_2( long New_Char,struct SARange & Start_Range, struct SARange & Dest_Range,BWT *fmi);
void Load_Indexes(BWT* & fwfmi, BWT* & revfmi,char *REVBWTINDEX,char *REVOCCFILE, char *REVSAFILE);

#endif
