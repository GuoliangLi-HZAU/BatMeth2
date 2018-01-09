#ifndef COMMON_H
#define COMMON_H
#define QLIMIT_FLOAT 30.0f
#define QLIMIT QLIMIT_FLOAT 
#define MULTIHIT_COUNT 3000
#define PLUS 50 
#define MINUS 51
#define SOLEXA_READS 0
#define SOLID_READS  1
#define SINGLE_END 1
#define PAIRED_END 2
#define MAX_MISMATCHES_BOUND 16
#define FQ	2
#define FA	3
#define MAXDES 500 
#define MAXTAG 500 
#define MAXBRANCH 500
#define REVERSE 1
#define FORWARD 0
#define SAINTERVAL 8
#define BRANCHTHRESHOLD 0//80 //30 //Threshold at which to check the BWT instead of branching
#define PAIR_END_SEPERATOR '\t'
#define PAIREND 2
#define INDELMARK (MAXTAG-1)
#define INSERTMARK 63
#define DELETEMARK 70
#define BATFILEMODE 	1	//pairing from already mapped batman file..
#define NORMALFILEMODE  2
#include <limits.h>
#include <float.h>
#include "zlib.h"
#include <map>
#include <vector> 
#include <queue>
#include <set>
#include <string>
#define NO 30
const int MAX_SIGLEN=50;
extern "C" 
{
	#include "iniparser-s.h"
	#include <time.h>
	#include "MemManager-s.h"
	#include "MiscUtilities-s.h"
	#include "TextConverter-s.h"
	#include "BWT-s.h"
}

enum {VERYFAST=0,FAST=1,SENSITIVE=2,VERYSENSITIVE=3};
enum {UNIQUEHIT/*only one hit*/,SHARP_UNIQUEHIT/*well resolved multi hit*/,MULTI_HIT/*Top hit, but multi hits*/,SW_RECOVERED/*same as above, but recovered with SW,no multihit indicated by Sub_Opt_Score=INT_MAX*/,SW_INDELSTEP/*multi hits found at indel stage*/,UNMAPPED/*No mapping*/,UNRESOLVED_HIT/*Several hits with the same score*/};
const int MX=6,MN=2,BOPEN=6,BEXT=3,MATCH_BONUS=0;//2;
//{-----------------------------  STRUCTS  -------------------------------------------------/

struct Ann_Info
{
	unsigned Size;
	unsigned Cumulative_Size;
	char* Name;
};

struct Offset_Record
{
	char Genome[40];
	unsigned Offset;
	int index;
};
struct Gene_Hash
{
	char* Genome;
	int Index;
};
struct OUTPUT//info for writing output..
{
	char SAM;
	char PLUSSTRAND;
	char MaxHits;
	char Offset;
	char* Buffer;
	char* Buffer_End;
	unsigned* Location_Array;
	Offset_Record* Genome_Offsets;
	int Genome_Count;
	int Length_Array[3];
	char FILETYPE;
};
/*
struct OUTPUT//info for writing output..
{
	char SAM;
	char PLUSSTRAND;
	char MaxHits;
	char Offset;
	char* Buffer;
	unsigned* Location_Array;
	Offset_Record* Genome_Offsets;
	int Genome_Count;
	int Length_Array[3];
	char FILETYPE;
};
*/
struct Header
{
	char ID[3] ;
	unsigned MAXHITS;
	char FILETYPE;
	char HITMODE;
	char IGNOREHEAD;
	char Index_Count;
	int Tag_Length;
	char Print_Desc;//print the tag desc?
}__attribute__((__packed__));

struct In_File
{
	int MAXHITS;
	int STRINGLENGTH;
	char FILETYPE;
	char MAX_MISMATCHES;
	int Stat_Size;
	int TAG_COPY_LEN;
	char ROLLOVER;
	char SCANBOTH;
	char PRINT_DESC;
	char LOADREVERSEONLY;
	char NORMAL_TAGS;
	int Length_Array[3];
	int HEAD_LENGTH,TAIL_LENGTH;
	char *Positive_Head,*Positive_Tail;
	char Tag_Copy[MAXTAG+1];
	off64_t File_Size;
};


// LEN stores the information related to stringlength...
struct LEN
{
	int STRINGLENGTH_ORG;//Real length of string
	int STRINGLENGTH;//effective length of string = Forced length-ignore head
	int STRINGLENGTHl;
	int STRINGLENGTHr; 
	int LH; 
	int RH; 
	int LHQL; 
	int LHQR; 
	int RHQL; 
	int RHQR; 
	int LHl; 
	int RHl; 
	int LHQLl; 
	int LHQRl; 
	int RHQLl; 
	int RHQRl; 
	int LHr; 
	int RHr; 
	int LHQLr; 
	int LHQRr; 
	int RHQLr; 
	int RHQRr;
	int IGNOREHEAD;
	bool batmeth;
};

struct READ
{
	char Description[MAXDES];
	char Raw_Tag_Copy[MAXDES];
	char Tag_Copy[MAXTAG];
	char Quality[MAXTAG];
	char Plus[MAXTAG];
	int NCount;//Number of N's
	char N[MAXTAG];
	char NLocations[MAXTAG];
	unsigned Tag_Number;//Head =1, Tail =2
	unsigned Read_Number;
	int Real_Len;
	off64_t FLength;
	char state;//0 raw state; 1 C2T; 2 G2A; 
};

struct FMFILES
{
	char* DISCORDANTFILE;
	char* PATTERNFILE; 
	char* PATTERNFILE1;
	char* HITSFILE;//=HITSFILE_DEF;//file to output hits...
	char* BWTFILE ; 
	char* OCCFILE ;
	char* REVBWTINDEX;
	char* REVOCCFILE;
	char* REVSAFILE;
	char* SAFILE;
	char* BINFILE;
	char* PACFILE;
	char* LOCATIONFILE ; 
	char* INPUTFILE;
	char* OUTPUTFILE;
	char* INDFILE;
	char* RANGEFILE;
	char* NLOCATIONFILE;
	char* BLKFILE;
	char* SINGLEFILE;

};

struct BATREAD
{
	int StringLength;
	int IGNOREHEAD;
	int NCount;
	char Forward[MAXTAG];
	char Complement[MAXTAG]; 
	char Forward_raw[MAXTAG];//moxian
	char Complement_raw[MAXTAG]; 
};

struct INFILE
{
	FILE* Input_File;
	off64_t File_Size;
	char FILETYPE;
	char SOLID;
	int TAG_COPY_LEN;
	char TAB_SEPERATED;
	char *Buffer;
	int PAIR_LENGTH_RIGHT;
	int PAIR_LENGTH_LEFT;
};

struct FILELIST
{
	FILE* Head;
	FILE* Tail;
};

struct BATPARAMETERS
{
	unsigned MAXHITS;
	int SCANMODE;
	char* PATTERNFILE;
	char* PATTERNFILE1;
	char* MISFILE1;
	char* MISFILE2;
	char* UNMAPPED_FILE;
	char* MULTI_FILE;
	char ONEFMINDEX;
	char UNMAPPED;
	int MAX_MISMATCHES;
	char PAIRING_MODE;
	unsigned MAX_TAGS_TO_PROCESS;
	int NTHREADS;
	int IGNOREHEAD; 
	int FORCELENGTH;
	int INSERTSIZE;
	int INDELSIZE;
	unsigned OFFSET;
	int PLUSSTRAND;
	int SW_FLANKSIZE;//Length of string to do S/W... A portion of this size will be taken from the either side of insert size + good hit.
	int FLANKSIZE;
	char LOADREVERSEONLY;
	char MISHITS;
	char USELOCATION; 
	int ORIENTATIONS;//how many strand combinations to try..
	int Patternfile_Count;//0=single end,1= paired end.. 
	int Misfile_Count;
	char ROLLOVER;
	char SCANBOTH;
	char SMITH_WATERMAN;
	char VERIFY;
	char LOG;
	char SOLIDMAP;
	int STD;
	char FORCESOLID;
	char* CMD_Buffer;
};

struct FastaSeq 
{
      int   len; /* the actual string length of seq */
      char* seq; /* the sequence itself */
};

typedef struct gz_stream {
    z_stream stream;
    int      z_err;   /* error code for last stream operation */
    int      z_eof;   /* set if end of input file */
    FILE     *file;   /* .gz file */
    Byte     *inbuf;  /* input buffer */
    Byte     *outbuf; /* output buffer */
    uLong    crc;     /* crc32 of uncompressed data */
    char     *msg;    /* error message */
    char     *path;   /* path name for debugging only */
    int      transparent; /* 1 if input file is not a .gz file */
    char     mode;    /* 'w' or 'r' */
    z_off_t  start;   /* start of compressed data in file (header skipped) */
    z_off_t  in;      /* bytes into deflate or inflate */
    z_off_t  out;     /* bytes out of deflate or inflate */
    int      back;    /* one character push-back */
    int      last;    /* true if push-back is last character */
} gz_stream;


struct HEADER
{
	char ID[3] ;
	unsigned MAXHITS;
	char FILETYPE;
	char HITMODE;
	char IGNOREHEAD;
	char Index_Count;
	int Tag_Length;
	char Print_Desc;//print the tag desc?
}__attribute__((__packed__));

struct SARange
{

	unsigned long Start;
	unsigned long End;
	int Level;//at what relative node are we?
	char Mismatches;//number of mismatches at this level
//	char Branches;//number of branches at this level
	unsigned char Skip;
	//int Tag;//Tag number
	char FMIndex;
	char Strand;
	unsigned Mismatch_Char;//2|2|...
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
//	unsigned char Branches_Pos[MAXBRANCH];//Branch position
        unsigned Mismatch_PosX;
//	unsigned long Start_Branch;
//	unsigned long End_Branch;
//	int DeleteStep;
//	int InsertStep;
	bool Indel;
//	int Level_Branch;
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
	unsigned Mismatch_PosX;//int?
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...
	unsigned Mismatch_Char;//2|2|...


}__attribute__((__packed__));

struct Mismatches_Record_GIS
{

	unsigned Mismatch_Char;//2|2|...
	unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];//BTS_PER_LOC|BTS_PER_LOC|...

}__attribute__((__packed__));


struct Output_Record
{
	unsigned Tag;
	unsigned  Start;
	char Index;
	unsigned char Skip;
	char Mismatches;
	int Gap;
}__attribute__((__packed__));

struct Branches
{
	long Is_Branch [4];
};

struct MEMLOOK
{
	int Lookupsize;
	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
};

struct MEMX
{
	SARange Branch_Ranges[4];

	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
	unsigned ARRAY_BOUND;
	unsigned END_BOUND;
	unsigned One_Mis_Bound;
	unsigned Hits;
	unsigned short Stats[7];

	char Guessed;
	char FMIndex;
	char Strand;
	char Stat_Size;
	char Larger_Than_Ten;
	char Extend;

	char* Write_Buffer;
	char* Current_Tag;
	char* Current_Tag_raw;

	int Hit_Array_Ptr;
	int Lookupsize;
	int Last_Mismatch_Written;
	SARange* Hit_Array;
	READ Read;
	OUTPUT Output;
	LEN L;

	int Left_Mishits_Pointer;
	int Right_Mishits_Pointer;
	int Possible_20_Pointer;
	int Possible_02_Pointer;
	int Mismatches_Forward_Pointer;
	int Mismatches_Backward_Pointer;
	int Two_Mismatches_At_End_Pointer;
	int Two_Mismatches_At_End_Forward_Pointer;
	int Possible_03_Pointer;
	int Possible_30_Pointer;
	int Possible_04_Pointer,Possible_40_Pointer,Possible_50_Pointer;
	int Possible_05_Pointer;
	int Mismatches_Forward_Pointer_Last4;
	int Left_Mishits_Pointer_1;
	int Mismatches_Forward_Pointer_Last5;

	int Least_Mis;
	int Best_Quality;

	SARange* BMHStack;
	SARange* FSHStack;
	SARange* FSHStackX0X;
	SARange* FSSStack;
	SARange* FSSStackX;
	SARange* BMStack;
	SARange* BMStackX;
	SARange* BMStack_X11;
	SARange* BMStack_X11H;
	SARange* PSBStack;

	SARange* Exact_Match_Forward;
	SARange* Exact_Match_Backward;
	SARange* Left_Mishits;
	SARange* Right_Mishits;
	SARange* Mismatches_Backward;
	SARange* Mismatches_Forward;
	SARange* Two_Mismatches_At_End_Forward;
	SARange* Two_Mismatches_At_End;
	SARange* Possible_20;
	SARange* Possible_02;
/*
	int Longest_H0_Pointer;
	int Longest_T0_Pointer;
*/
	int Longest_H1_Pointer;
	int Longest_T1_Pointer;
	
	SARange* Longest_H0;
	SARange* Longest_T0;
	SARange* Longest_H1;
	SARange* Longest_T1;

};

struct GUESS
{
	MEMX Guessed;
	char* Guessed_Read;
	MEMX Guess_Complement;
	char* Guessed_ReadC;
};
struct MEM
{
	unsigned* Forward_Start_LookupX;
	unsigned* Forward_End_LookupX;
	unsigned* Backward_Start_LookupX;
	unsigned* Backward_End_LookupX;
	unsigned ARRAY_BOUND;
	unsigned END_BOUND;

	char* Write_Buffer;
	int Lookupsize;

	SARange* BMHStack;
	SARange* FSHStack;
	SARange* FSHStackX0X;
	SARange* FSSStack;
	SARange* FSSStackX;
	SARange* BMStack;
	SARange* BMStackX;
	SARange* BMStack_X11;
	SARange* BMStack_X11H;
	SARange* PSBStack;

	SARange* Exact_Match_ForwardF;
	SARange* Exact_Match_BackwardF;
	SARange* Left_MishitsF;
	SARange* Right_MishitsF;
	SARange* Mismatches_BackwardF;
	SARange* Mismatches_ForwardF;
	SARange* Two_Mismatches_At_End_ForwardF;
	SARange* Two_Mismatches_At_EndF;
	SARange* Possible_20F;
	SARange* Possible_02F;

	SARange* Exact_Match_ForwardC;
	SARange* Exact_Match_BackwardC;
	SARange* Left_MishitsC;
	SARange* Right_MishitsC;
	SARange* Mismatches_BackwardC;
	SARange* Mismatches_ForwardC;
	SARange* Two_Mismatches_At_End_ForwardC;
	SARange* Two_Mismatches_At_EndC;
	SARange* Possible_20C;
	SARange* Possible_02C;

};



struct Thread_Arg
{
	BWT* fwfmi;
	BWT* revfmi;
	FILE* Output_File;
	LEN  L; 
	char ONEFMINDEX;
	MEMLOOK MLook;
	OUTPUT Output;
	BATPARAMETERS Bat_Para;
	INFILE In;
	int ThreadID;
};


struct Threading
{
	pthread_t Thread;
	unsigned r;
	void *ret;
	Thread_Arg Arg;
};

struct Alignment
{
	unsigned Loc;
	char Sign;
	char *Chr;
	char source;
	//unsigned char Bad_Loc[SW_MAX_MIS_TO_STORE];
	//char Bad_Char[10];
	//char Edit_Distance;
	int QualityScore;
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
	int Top_Penalty;
        //
	int Clip;
	char Extend;
	bool Rescued;
	bool Do_Rescue;
};

class Comp_Alignment
{
	public:
		bool operator()(Alignment& A1,Alignment& A2)
		{
			if (A1.Score <= A2.Score) return true; else return false;
		}
};
struct Align_Hit
{
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> H0; 
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> H1; 
        std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> H2; 
        
	std::priority_queue <Alignment,std::vector <Alignment>,Comp_Alignment> AlignTmp;
        std::set<unsigned> AM;//Align_Maps
};
struct Hit_Info
{
	char* MH;
	char* Chr;
	char Status;
	unsigned Loc;
	unsigned Org_Loc;
	unsigned Sub_Opt_Hit;
	char Sign;
	int Score;
	int Sub_Opt_Score;
	int QScore;
	int SW_Score;
	int SW_Sub_Opt_Score;
	int Hits;
	int Mismatch;
	int Indel;
	float BQScore;
	float Sub_BQScore;
	char Cigar[MAX_SIGLEN+1];
	int Clip_H,Clip_T;
};

struct Hit_
{
	char* Chr;
	unsigned Loc;
	char Sign;
};


struct Hit_Status//info about hits found in a scanning phase for head or tail..
{
	int Last_Mis;
	SARange* Sub_Opt_StartF;//start of suboptimal hits..
	SARange* Sub_Opt_StartC;
	int Tophits;
	int Subopthits;
};

//jq
struct Scan_Params
{
	int Num_MM;
	int Num_GAP;
	int Cost;
};
//}-----------------------------  STRUCTS  -------------------------------------------------/*
#include <signal.h>
#define Breakpoint raise(SIGINT)

#endif
