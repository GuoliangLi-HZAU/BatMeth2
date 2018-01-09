#include "common-s.h"
#include "rqindex-s.h"

#ifdef __MAIN_CODE__
	#ifndef __MAIN_INSERTED__
	#define __MAIN_INSERTED__
		#ifdef __cplusplus
		std::map <unsigned, Ann_Info> Annotations;
		std::map <unsigned, Ann_Info> ::iterator S,E;
		#endif 
		bool non_directional=false;
		unsigned MAX_Hits=1000;
		bool FindInDels=true;
		bool List_Filter = true;
		BWT *fwfmiCT,*revfmiCT,*fwfmiGA,*revfmiGA;//moxian
		//BWT *fwfmi,*revfmi; //moxian
		RQINDEX RQ,RQ_CT,RQ_GA/*,RQ_revGA*/;
		RQINDEX RQHALF,RQHALF_CT,RQHALF_GA;
		unsigned Entries,Entries_CT,Entries_GA,Entries_revGA;//
		unsigned Entries_Half,Entries_Half_CT,Entries_Half_GA;
		//unsigned Conversion_Factor;
		//int Genome_Count=0;
		char* CONFIG_FILE=NULL;
		char* LOG_SUCCESS_FILE=NULL;
		char* REJECT_FILE=NULL;
		char Multi_HitH[MULTIHIT_COUNT];
		char Multi_HitT[MULTIHIT_COUNT];
		//char* Temp_Disc_Hit_Buffer;
		//char* Disc_Hit_Buffer_Org;

		FILE* Mishit_File1;
		FILE* Mishit_File2;
		FILE* Log_SFile;
		FILE* Reject_File;
		FILE* GLog_File;
		bool Only_non_directional=false;
		unsigned char* Original_Text_Ori;
		unsigned char* Original_Text_CT;
		unsigned char* Original_Text_GA;
		//unsigned char* Original_Text_revGA;
		const char* Code_To_Char="acgt";
		const char* Code_To_CharU="ACGTN";
		char Code_To_CodeC[5];

		//Offset_Record Genome_Offsets[80];
		//unsigned Location_Array[80];
		int LEAST_MIS=6;

		float SW_EDIT_PERC=0;
		int SW_EDIT_DIST_MAX=9;
		int TRIM_LENGTH=0;
		int NCOUNT=3;
		int Lens=75;
		int GOOD_N_LEN=20;
		int ACC_SCORE=1;

	//mismatch_max
		int MAX_MISMATCHES=5;
		char SWBUG=FALSE;
		char UNIQ_HIT=FALSE;
		char SW_EDIT_DIST_MAX_CHANGED=FALSE; 
		char SW_IS_LOW=FALSE;
		char PRINT_ONLY_UNIQUE_DISC=FALSE;
		char CLEAN_HITS_ONLY=FALSE;
		char DICT=FALSE;
		char PROGRESSBAR=TRUE;
		char MISC_VERB=TRUE;
		char RECOVER_N=TRUE;
		unsigned SW_HITS,SW_DUP;
		unsigned Total_Paired=0;
		unsigned Total_Mapped=0;
		unsigned Total_Reads=0;
		int READS_TO_PROCESS =0;
		int PAIRING_SCHEME=5335;
		bool PRINT_UNIQUE_SW_ONLY=0;//true; 
		bool PRINT_TOP_HIT=1;//false;
		bool SUPOPT_EXIST_FILTER=0;
		bool PHRED_FILTER=true;
		int HITS_IN_SAM=INT_MAX;//40;
		bool HEURISTIC=0;
		bool SKIPHARD=false;//skips time consuming middle indel scan..
		bool PRINT_MISHIT=false;
		bool LAZY_PHRED=false;
		bool BAYES_PHRED=true;
		bool PHRED64=false;

		//char Description[MAXDES+1];
		//char Default_Cigar[50];
		char Char_To_Code[255];
		char Char_To_CodeC[256];
		char Char_To_CharC[256];

		//char Asterisk[]="*";
		//char Empty_String[]="";
		char GENOMEFILE_DEFAULT[]="genome";
		char PATTERNFILE_DEFAULT[]="tags.fq";
		char HITSFILE_DEFAULT[]="hits.txt";
		char UNIQUEFILE_DEFAULT[]="unique.txt";
		char MISHITFILE_DEFAULT[]="mishits.fq";
		char BWTFILE_DEFAULT[] = "genome.bwt";//"genome.bwt";// 
		char OCCFILE_DEFAULT[] ="genome.fmv";//"genome.fmv";//
		char REVBWTINDEX_DEFAULT[] ="revgenome.bwt";//"revgenome.bwt";//
		char REVOCCFILE_DEFAULT[] ="revgenome.fmv";//"revgenome.fmv";//
		char REVSAFILE_DEFAULT[] ="revgenome.sa";
		char SAFILE_DEFAULT[] ="genome.sa";
		char BINFILE_DEFAULT[] = "genome.bin";//"genome.bwt";// 
		char INPUTFILE_DEFAULT[]="hits.txt";
		char OUTPUTFILE_DEFAULT[]="output.txt";
		char LOGFILE[]="run.log";
		char LOCATIONFILE_DEFAULT[]="location";
		char OUTFILE[]="-";//"pairings.txt";
		char HITSFILE_DEF[]="hits.txt";

		char SOLID;
		char RESCUE_DISC=FALSE;
		char BWAPOL=FALSE;
		char NPOLICY=0;
		char WRITE_DISCORDANT=FALSE;
		char WRITE_SINGLE=FALSE;
		char NFILE=FALSE;
		char BEST=FALSE;
		char WRITE_UNMAPPED=FALSE;
		char WRITE_MULTI=FALSE;
		char PRIORITYMODE=FALSE;
		char PRINT_DICTIONARY=TRUE;
		char INPUT_FILE_TYPE;
		float Mis_FA_Score,Match_FA_Score; 

		//jq
		int QP[101];

		int RECLIMIT;
		int QUAL_START;
		int SW_MIN_MATCH_LEN =20;//6;//20
		int Missed_Hits=0;
		int MAX_SW=INT_MAX;//5000;//120;//5000;//120//1;//MAXCOUNT_DEFAULT; 
		#ifndef NDEBUG
		char* DISC_HIT_BUFFER_ORG;
		#endif
		int32_t match, mismatch, gap_open, gap_extension, path, n, filter;
		//int8_t* mata;
	#endif
#else 
#ifndef __GLOBAL_VAR__
#define __GLOBAL_VAR__
#ifdef __cplusplus
extern std::map <unsigned, Ann_Info> Annotations;
extern std::map <unsigned, Ann_Info> ::iterator S,E;
#endif 
extern int MAX_MISMATCHES;
extern BWT *fwfmiCT,*revfmiCT,*fwfmiGA,*revfmiGA;//moxian
extern RQINDEX RQ,RQ_CT,RQ_GA/*,RQ_revGA*/;
extern RQINDEX RQHALF,RQHALF_CT,RQHALF_GA;
extern unsigned Entries,Entries_CT,Entries_GA/*,Entries_revGA*/;//
extern unsigned Entries_Half,Entries_Half_CT,Entries_Half_GA;
//extern unsigned Conversion_Factor;
//extern int Genome_Count;
extern char* CONFIG_FILE;
extern char* LOG_SUCCESS_FILE;
extern char* REJECT_FILE;
extern char Multi_HitH[MULTIHIT_COUNT];
extern char Multi_HitT[MULTIHIT_COUNT];
//extern char* Temp_Disc_Hit_Buffer;
//extern char* Disc_Hit_Buffer_Org;

extern FILE* Mishit_File1;
extern FILE* Mishit_File2;
extern FILE* Log_SFile;
extern FILE* Reject_File;
extern FILE* GLog_File;

extern unsigned char* Original_Text_Ori;
extern unsigned char* Original_Text_CT;
extern unsigned char* Original_Text_GA;
//extern unsigned char* Original_Text_revGA;
extern const char* Code_To_Char;
extern const char* Code_To_CharU;
extern char Code_To_CodeC[5];
extern LEN GlobL;

//extern Offset_Record Genome_Offsets[80];
//extern unsigned Location_Array[80];
extern int LEAST_MIS;

extern float SW_EDIT_PERC;
extern int SW_EDIT_DIST_MAX;
extern int TRIM_LENGTH;
extern int NCOUNT;
extern int Lens;
extern int GOOD_N_LEN;
extern int ACC_SCORE;

extern char SWBUG;
extern char UNIQ_HIT;
extern char SW_EDIT_DIST_MAX_CHANGED; 
extern char SW_IS_LOW;
extern char PRINT_ONLY_UNIQUE_DISC;
extern char CLEAN_HITS_ONLY;
extern char DICT;
extern char PROGRESSBAR;
extern char MISC_VERB;
extern char RECOVER_N;
extern unsigned SW_HITS,SW_DUP;
extern unsigned Total_Paired;
extern unsigned Total_Mapped;
extern unsigned Total_Reads;
extern int READS_TO_PROCESS ;
extern int PAIRING_SCHEME;
extern bool PRINT_UNIQUE_SW_ONLY;//true; 
extern bool PRINT_TOP_HIT;//false;
extern bool SUPOPT_EXIST_FILTER;
extern bool PHRED_FILTER;
extern int HITS_IN_SAM;
extern bool HEURISTIC;
extern bool SKIPHARD;//skips time consuming middle indel scan..
extern bool PRINT_MISHIT;
extern bool LAZY_PHRED;
extern bool BAYES_PHRED;
extern bool PHRED64;

//char Description[MAXDES+1];
//extern char Default_Cigar[50];
extern char Char_To_Code[255];
extern char Char_To_CodeC[256];
extern char Char_To_CharC[256];

//extern char Asterisk[];
//extern char Empty_String[];
extern char GENOMEFILE_DEFAULT[];
extern char PATTERNFILE_DEFAULT[];
extern char HITSFILE_DEFAULT[];
extern char UNIQUEFILE_DEFAULT[];
extern char MISHITFILE_DEFAULT[];
extern char BWTFILE_DEFAULT[];//"genome.bwt";// 
extern char OCCFILE_DEFAULT[];//"genome.fmv";//
extern char REVBWTINDEX_DEFAULT[];//"revgenome.bwt";//
extern char REVOCCFILE_DEFAULT[];//"revgenome.fmv";//
extern char REVSAFILE_DEFAULT[];
extern char SAFILE_DEFAULT[];
extern char BINFILE_DEFAULT[];//"genome.bwt";// 
extern char INPUTFILE_DEFAULT[];
extern char OUTPUTFILE_DEFAULT[];
extern char LOGFILE[];
extern char LOCATIONFILE_DEFAULT[];
extern char OUTFILE[];
extern char HITSFILE_DEF[];

extern char SOLID;
extern char RESCUE_DISC;
extern char BWAPOL;
extern char NPOLICY;
extern char WRITE_DISCORDANT;
extern char WRITE_SINGLE;
extern char NFILE;
extern char BEST;
extern char WRITE_UNMAPPED;
extern char WRITE_MULTI;
extern char PRIORITYMODE;
extern char PRINT_DICTIONARY;
extern char INPUT_FILE_TYPE;
extern float Mis_FA_Score,Match_FA_Score; 

//jq
extern int QP[101];

extern int RECLIMIT;
extern int QUAL_START;
extern int SW_MIN_MATCH_LEN;//6;//20
extern int Missed_Hits;
extern int  MAX_SW;
#ifndef NDEBUG
extern char* DISC_HIT_BUFFER_ORG;
#endif
extern int32_t match, mismatch, gap_open, gap_extension, path, n, filter;
//extern int8_t* mata;
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------
#endif
#endif
