#include <getopt.h> 
#include "string.h" 
#include "command.h"
#include "common.h"
#include "global.h"
//{---------------------------- Command Line  -------------------------------------------------
extern int DEB;
extern int LENGTH_CUTOFF;
extern int SW_SIMILARITY_FOR_RESCUE;
extern int INSERTSIZE;
extern int STD;
extern int Inter_MM;
extern int Max_MM;
extern int Max_MM_GAP;
extern int QUALITYCONVERSIONFACTOR;
extern int THREAD;
extern bool FASTSW;
extern bool USE_MULTI_OUT;
extern bool STACK_LOWQ;
extern bool FASTDECODE;
extern int BOOST;
extern int TOP_TEN;
extern int MAX_PER_LIST;
extern int SCOREGAP;
extern bool DEBUG_SEGS;
extern bool PAIRED;
extern bool Hard_Penalty;
extern bool REALN;
extern bool non_directional;
extern int SMALLSEED;
//extern bool TESTMODE;

option Long_Options[]=
{
{"print_unique_sw_only",0,NULL,1},
{"print_tophits",0,NULL,2},
{"heuristic",0,NULL,3},
{"no_phred_filter",0,NULL,4},
{"easyphred",0,NULL,5},
{"matchscore",1,NULL,6},
{"mismatchscore",1,NULL,7},
{"gapopen",1,NULL,8},
{"gapextend",1,NULL,9},
{"general",0,NULL,10},
{"swlimit",1,NULL,11},
{"intermm",1,NULL,12},
{"maxmm",1,NULL,13},
{"maxmmgap",1,NULL,14},
{"illumina",0,NULL,15},
{"threads",1,NULL,'p'},
{"splitout",0,NULL,17},
{"stacklowq",0,NULL,18},
{"indelsize",1,NULL,19},
{"mode",1,NULL,20},
{"fastdecode",0,NULL,21},
{"boost",1,NULL,22},
{"nofastsw",0,NULL,23},
{"tophits",1,NULL,24},
{"listdepth",1,NULL,25},
{"scoregap",1,NULL,26},
{"debug",0,NULL,27},
{"concthreshold",1,NULL,28},
{"lengthcutoff",1,NULL,29},
{"softpenalty",0,NULL,30},
{"norealign",0,NULL,31},
{"seed",1,NULL,32},
{"help",0,NULL,'h'},
{"output",1,NULL,'o'},
{"genome",1,NULL,'g'},
{"verbose",1,NULL,'v'},
{"noheader",0,NULL,'H'},
{"buffersize",1,NULL,'b'},
{"insertsize",1,NULL,'s'},
{"maxmismatches",1,NULL,'n'},
{"editdist",1,NULL,'e'},
{"matchlen",1,NULL,'E'},
{"bwapol",0,NULL,'B'},
{"flanksize",1,NULL,'f'},
{"maxtags",1,NULL,'t'},
{"swon",0,NULL,'w'},
{"swlow",0,NULL,'Z'},
{"ntorand",0,NULL,'N'},
{"ncount",1,NULL,'C'},
{"outputfile",1,NULL,'o'}, 
{"maxhits",1,NULL,'m'},
{"location",optional_argument,NULL,'l'},
{"rollover",0,NULL,'r'},
{"std",1,NULL,'d'},
{"unmapped",1,NULL,'X'},
{"single",1,NULL,'u'},
{"disc",1,NULL,'D'},
{"nhits",1,NULL,'r'},
{"multihit",1,NULL,'M'},
{"forcelength",1,NULL,'F'},
{"unique",1,NULL,'U'},
{"lowestmis",1,NULL,'P'},
{"logfile",1,NULL,'L'},
{"non_directional",0,NULL,16},
{0,0,0,0}
};
//}---------------------------- Command Line -------------------------------------------------
//{-----------------------------  Parse Command Line  -------------------------------------------------
void usage()
{
        fprintf(stderr,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
        fprintf(stderr, "Program: BatMeth2 \n");
        fprintf(stderr, "Version: v1.0\n\n");
        fprintf(stderr,"[ Single-end-reads ]\n");
        fprintf(stderr,"Usage:     batmeth2 -g INDEX -i INPUT -o OUTPUT\n");
        fprintf(stderr,"Example:   batmeth2 -g /data/index/hg19/hg19.fa -i Read.fq -o outPrefix -p 6 -n 2\n\n");
        fprintf(stderr,"[ Paired-end-reads ]\n");
        fprintf(stderr,"Usage:     batmeth2 -g INDEX -i INPUT_left -i INPUT_right -o OUTPUT\n");
        fprintf(stderr,"Example:   batmeth2 -g /data/index/hg19/hg19.fa -i Read_R1_left.fq -i Read_R2_right.fq -o outPrefix -p 6 -n 2\n\n");
//      fprintf(stderr,"----------------------------------------------------------------------\n\n");
        fprintf(stderr,"Parameters : \n");
        fprintf(stderr,"      --inputfile | -i <filename>   Name of input file\n");
        fprintf(stderr,"      --genome | -g <filename>      Name of the genome mapped against\n");
        fprintf(stderr,"      --outputfile | -o <filename>  Name of output file prefix\n");
        fprintf(stderr,"      --buffersize | -b <integer>   Size of disk buffers\n");
        fprintf(stderr,"      --indelsize                   indel size\n");
	fprintf(stderr,"      --seed                        seed length for alignment. suggest: n/2. n = length of the read.\n");
        fprintf(stderr,"      --maxhits | -m                maximum number of pairings to find\n");
        fprintf(stderr,"      --maxmismatches | -n          maximum mismatches allowed due to seq. errors\n");
        fprintf(stderr,"      --non_directional             Alignments to all four bisulfite strands will be reported. Default: OFF. \n");
        fprintf(stderr,"      --single | -u                 output singly mapped hits to file \n");
        fprintf(stderr,"      --editdist | -e               maximum edit distance/percentage  to allow in mate rescue\n");
        fprintf(stderr,"      --matchlen | -E               minimum match len to allow in mate rescue\n");
        fprintf(stderr,"      --maxtags | -t <integer>      maximum pairs to scan\n");
        fprintf(stderr,"      --std | -d <integer>          standard deviatiion of reads distribution\n");
        fprintf(stderr,"      --swon | -w(=rescdisc)        perform Smith-Waterman mate rescue. if rescdisc is set, try to rescue discordant pairs before marking them as discordant\n");
        fprintf(stderr,"      --swlow | -Z                  Assign low quality to rescued mate\n");
        fprintf(stderr,"      --flanksize | -f <integer>    size of flanking region for Smith-Waterman\n");
        fprintf(stderr,"      --forcelength | -F <integer>  map only the first <integer> bases\n");
        fprintf(stderr,"      --bwapol | -B                 descend only one mismatch mote than opt hit\n");
        fprintf(stderr,"      --unique | -U                 Mark unique hits\n");
        fprintf(stderr,"      --lowestmis | -P              Pair Hits having lowest count of mismatches in Head and Tail\n");
        fprintf(stderr,"      --config | -c <file>          Load settings from configuration <file>  \n");
        fprintf(stderr,"      --verbose | -v <number>       1-turn off general info 2-turn off progress bar \n");
        fprintf(stderr,"      --ntorand | N                 Replace N's in the read with random characters\n");
        fprintf(stderr,"      --ncount | C                  Maximum N's allowed in the read\n");
        fprintf(stderr,"      --logfile | -L                Log exit status to this log \n");
        fprintf(stderr,"      --heuristic                   Skip 2 mismatch scan of end anchors\n");
        fprintf(stderr,"      --no_phred_filter             do not apply phred filter\n");
        fprintf(stderr,"      --easyphred                   apply lenient phred filter\n");
        fprintf(stderr,"      --general                     try to get most accurate mapping\n");
        fprintf(stderr,"      --swlimit | <integer>         try at most <integer> sw extensions\n");
        fprintf(stderr,"      --threads | -p <integer>      Launch <integer> threads\n");
        fprintf(stderr,"      --NoInDels | -I               not to find the indels result\n");
        fprintf(stderr,"      --help | -h                   Print help\n\n");
        fprintf(stderr,"Note: To use BatMeth2, you need to first index the genome with `build_all genome'.\n\n");
}
void Parse_Command_line(int argc, char* argv[],char* & GENOME,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP)
{
	int Current_Option=0;
	const char* Short_Options ="C:c:e:E:d:D:f:F:x:X:s:ShHi:b:l::L:o:g:Pq:r:t:m:M:n:NbBw::R:Uu:v:Zj:p:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	BP.INDELSIZE=8;/*8|21*/BP.SCANMODE=SENSITIVE;MAX_SW=INT_MAX;
	match = 1;//1;
	mismatch = 2; //4;//4;
	gap_open = 6; //6;//11;//3;
	gap_extension = 1;
	char* Name;int Last_Dash;char* Genome_Name;
	char* Command_Line_Buffer;
	GENOME=0;

	Command_Line_Buffer=BP.CMD_Buffer;
	for(int i = 0; i < argc; i++) Command_Line_Buffer+=sprintf(Command_Line_Buffer," %s", argv[i]);sprintf(Command_Line_Buffer,"\n");

	Command_Line_Buffer=(char*)malloc(6000);
	int Par_Count=0,Verb,ti;
	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		Par_Count++;
		switch(Current_Option)
		{
			case 1:
				PRINT_UNIQUE_SW_ONLY=1;
				break;
			case 2:
				PRINT_TOP_HIT=1;
				break;
			case 3://heuristic
				BEST=FALSE;
				PHRED_FILTER=true;
				LAZY_PHRED=false;
				BAYES_PHRED=true;
				Max_MM_GAP=0;
				MAX_SW=100;
				HEURISTIC=1;
				break;
			case 4://no_phred_filter
				PHRED_FILTER=false;
				LAZY_PHRED=false;
				BAYES_PHRED=true;
				MAX_SW=5000;
				break;
			case 5:
				PHRED_FILTER=true;
				LAZY_PHRED=true;
				BAYES_PHRED=false;
				break;
			case 6:
				match=atoi(optarg);
				break;
			case 7:
				mismatch=atoi(optarg);
				break;
			case 8:
				gap_open=atoi(optarg);
				break;
			case 9:
				gap_extension=atoi(optarg);
				break;
			case 10://general
				SKIPHARD=true;
				break;
			case 11://swlimit
				if(!strcmp(optarg,"all"))
				{
					printf("FULL SW...\n");
					MAX_SW=INT_MAX;
				}
				else
					MAX_SW=atoi(optarg);
				break;
			case 12://intermm
				Inter_MM=atoi(optarg);
				break;
			case 13://maxmm
				Max_MM=atoi(optarg);
				break;
			case 14://maxmmgap
				Max_MM_GAP=atoi(optarg);
				break;
			case 15://illumina
				PHRED64=true;
				QUALITYCONVERSIONFACTOR=64;
				break;
			case 'p'://thread
				THREAD=atoi(optarg);
				break;
			case 16:
				non_directional=true;
				break;
			case 17://splitout split threaded out to many files..
				USE_MULTI_OUT=true;
				break;
			case 18://stacklowq 
				STACK_LOWQ=true;
				break;
			case 19://indelsize 
				BP.INDELSIZE=atoi(optarg);
				break;
			case 20://mode 
				if(!strcmp(optarg,"sensitive"))
				{
					BP.SCANMODE=SENSITIVE;
				}
				else if(!strcmp(optarg,"vsensitive"))
				{
					BP.SCANMODE=VERYSENSITIVE;
				}
				else if(!strcmp(optarg,"fast"))
				{
					BP.SCANMODE=FAST;
				}
				else if(!strcmp(optarg,"vfast"))
				{
					BP.SCANMODE=VERYFAST;
				}
				break;
			case 21://fastdecode.. 
				FASTDECODE=true;
				break;
			case 22://boost.. 
				BOOST=atoi(optarg);
				printf("BOOST Level :%d",BOOST);
				break;
			case 23://nofastsw.. 
				FASTSW=false;
				printf("Fast Smith-Waterman OFF\n");
				break;
			case 24://tophits.. 
				TOP_TEN=atoi(optarg);
				break;
			case 25://listdepth
				if(!strcmp(optarg,"all"))
				{
					printf("FULL DEPTH...\n");
					MAX_PER_LIST=INT_MAX;
				}
				else
					MAX_PER_LIST=atoi(optarg);
				break;
			case 26://scoregap
				SCOREGAP=atoi(optarg);
				break;
			case 27://debugsegs
				DEBUG_SEGS=true;
				DEB=true;
				break;
			case 28://concthreshold
				SW_SIMILARITY_FOR_RESCUE=atoi(optarg);
				break;
			case 29://lengthcutoff
				LENGTH_CUTOFF=atoi(optarg);
				break;
			case 30://softpenalty
				Hard_Penalty=false;
				break;
			case 31://norealign
				REALN=false;
				break;
			case 32: //seed
				SMALLSEED=atoi(optarg);
			case 'j':
				REJECT_FILE=optarg;
				Reject_File=File_Open(REJECT_FILE,"w");
				break;
			case 'm':
				if(!strcmp(optarg,"all")) MAXCOUNT=INT_MAX;
				else MAXCOUNT=atoi(optarg);
				break;
			case 'v':
				Verb=atoi(optarg);
				if (Verb == 1){MISC_VERB=FALSE;} else if (Verb ==2){PROGRESSBAR=0;}
				break;
			case 'F':
				TRIM_LENGTH=atoi(optarg);
				break;
			case 'Z':
				SW_IS_LOW=TRUE;
				break;
			case 'C':
				NCOUNT=atoi(optarg);
				break;
			case 'N':
				NPOLICY=1;
				break;
			case 'H':
				PRINT_DICTIONARY=FALSE;
				break;
			case 't':
				READS_TO_PROCESS=atoi(optarg);
				break;
			case 'r':
				if(!strcmp(optarg,"all")) HITS_IN_SAM=INT_MAX;
				else if((HITS_IN_SAM=atoi(optarg))<=0) {printf ("Parse_Command_line():Invalid parameter for -r\n");exit(-1);}
				break;
			case 'U':
				UNIQ_HIT=TRUE;
				if(MAXCOUNT==1) MAXCOUNT =2;
				break;
			case 'x':
				PRINT_MISHIT=true;
				if(optarg)
				{
					BP.MISFILE1=optarg;
				}
				else
				{
					BP.MISFILE1="mishit.fq";
				}
				break;
			case 'X':
				WRITE_UNMAPPED=TRUE;
				BP.UNMAPPED_FILE=optarg;
				break;
			case 'c':
				CONFIG_FILE=optarg;
				break;
			case 'M':
				WRITE_MULTI=TRUE;
				BP.MULTI_FILE=optarg;
				break;
			case 'e':
				if (atoi(optarg)==atof(optarg)) {SW_EDIT_DIST_MAX =atoi(optarg);SW_EDIT_DIST_MAX_CHANGED=TRUE;} 
				else SW_EDIT_PERC=atof(optarg);
				break;
			case 'E':
				SW_MIN_MATCH_LEN = atoi(optarg);
				break;
			case 'R':
				RECLIMIT = atoi(optarg);
				break;
			case 'P':
				PRIORITYMODE=TRUE;
				break;
			case 'B':
				BWAPOL=TRUE;
				break;
			case 'n':
				BP.MAX_MISMATCHES = atoi(optarg);
				Inter_MM=atoi(optarg);
				MAX_MISMATCHES=atoi(optarg);
				if (BP.MAX_MISMATCHES <0 or BP.MAX_MISMATCHES >MAX_MISMATCHES_BOUND) BP.MAX_MISMATCHES=2;
				break;
			case 'w':
				BP.SMITH_WATERMAN=TRUE;
				if (optarg && !strcmp(optarg,"rescdisc")) RESCUE_DISC=TRUE;
				break;
			/*case 'X':
				BP.ORIENTATIONS=4;
				BP.SCANBOTH=TRUE;
				break;*/
			case 'd':
				STD=atoi(optarg);
				break;
			case 'f':
				BP.SW_FLANKSIZE=atoi(optarg);
				break;
			case 'L':
				LOG_SUCCESS_FILE=optarg;
				Log_SFile=File_Open(LOG_SUCCESS_FILE,"w");
				break;
			case 'i':
				if(!BP.Patternfile_Count){BP.PATTERNFILE=optarg;}
				else 
				{
					BP.PATTERNFILE1=optarg;
					PAIRED=TRUE;
				}
				BP.Patternfile_Count++;
				BP.PAIRING_MODE=NORMALFILEMODE;
				break;
			case 'h':
				usage();exit(0);
			case 'q':
				BP.PAIRING_MODE=BATFILEMODE;
				BP.PATTERNFILE=optarg;
				break;
			case 'o':
				F.OUTPUTFILE=optarg;
				break;
			case 'D':
				F.DISCORDANTFILE=optarg;
				WRITE_DISCORDANT=TRUE;
				break;
			case 'u':
				F.SINGLEFILE=optarg;
				WRITE_SINGLE=TRUE;
				break;
			case 's':
				INSERTSIZE=atoi(optarg);
				break;
			case 'S':
				BP.FORCESOLID=TRUE;
				break;
			case 'l':
				if (optarg)
				{
					if (!strcmp(optarg,"N")) NFILE=TRUE;
					else F.LOCATIONFILE=optarg;	
				}
				BP.USELOCATION=TRUE;
				break;
			case 'g':
				Name=optarg;GENOME=optarg;Genome_Name=optarg;
				Last_Dash=0;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}
/*
				F.REVBWTINDEX = (char*)Command_Line_Buffer;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(F.REVBWTINDEX,optarg,Last_Dash);
				F.REVBWTINDEX[Last_Dash+0]='r';F.REVBWTINDEX[Last_Dash+1]='e';F.REVBWTINDEX[Last_Dash+2]='v';
				strcpy(F.REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(F.REVBWTINDEX+Last_Dash+3,".bwt"); 

				F.BWTFILE=F.REVBWTINDEX+500;
				strncpy(F.BWTFILE,optarg,Last_Dash);
				strcpy(F.BWTFILE+Last_Dash,Genome_Name+1);
				strcat(F.BWTFILE+Last_Dash,".bwt"); 


				F.REVOCCFILE = F.BWTFILE+500;
				strncpy(F.REVOCCFILE,optarg,Last_Dash);
				F.REVOCCFILE[Last_Dash+0]='r';F.REVOCCFILE[Last_Dash+1]='e';F.REVOCCFILE[Last_Dash+2]='v';
				strcpy(F.REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(F.REVOCCFILE+Last_Dash+3,".fmv"); 


				F.OCCFILE=F.REVOCCFILE+500;			
				strncpy(F.OCCFILE,optarg,Last_Dash);
				strcpy(F.OCCFILE+Last_Dash,Genome_Name+1);
				strcat(F.OCCFILE+Last_Dash,".fmv"); 

				F.SAFILE=F.OCCFILE+500;			
				strncpy(F.SAFILE,optarg,Last_Dash);
				strcpy(F.SAFILE+Last_Dash,Genome_Name+1);
				strcat(F.SAFILE+Last_Dash,".sa");

				F.REVSAFILE = F.SAFILE+500;
				strncpy(F.REVSAFILE,optarg,Last_Dash);
				F.REVSAFILE[Last_Dash+0]='r';F.REVSAFILE[Last_Dash+1]='e';F.REVSAFILE[Last_Dash+2]='v';
				strcpy(F.REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(F.REVSAFILE+Last_Dash+3,".sa"); 

				F.BINFILE=F.REVSAFILE+500;			
				strncpy(F.BINFILE,optarg,Last_Dash);
				strcpy(F.BINFILE+Last_Dash,Genome_Name+1);
				strcat(F.BINFILE+Last_Dash,".pac");

				if(!BP.USELOCATION)
				{

					F.LOCATIONFILE=F.BINFILE+500;			
					strncpy(F.LOCATIONFILE,optarg,Last_Dash);
					strcpy(F.LOCATIONFILE+Last_Dash,Genome_Name+1);
					strcat(F.LOCATIONFILE+Last_Dash,".ann.location");
				}

                                F.BLKFILE = F.LOCATIONFILE+500;
                                strncpy(F.BLKFILE,optarg,Last_Dash);
                                strcpy(F.BLKFILE+Last_Dash,Genome_Name+1);
                                strcat(F.BLKFILE+Last_Dash,".blk.");

                                F.INDFILE = F.BLKFILE+500;
                                strncpy(F.INDFILE,optarg,Last_Dash);
                                strcpy(F.INDFILE+Last_Dash,Genome_Name+1);
                                strcat(F.INDFILE+Last_Dash,".ind.");

                                F.RANGEFILE = F.INDFILE+500;
                                strncpy(F.RANGEFILE,optarg,Last_Dash);
                                strcpy(F.RANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(F.RANGEFILE+Last_Dash,".range");

                                F.NLOCATIONFILE = F.RANGEFILE+500;
                                strncpy(F.NLOCATIONFILE,optarg,Last_Dash);
                                strcpy(F.NLOCATIONFILE+Last_Dash,Genome_Name+1);
                                strcat(F.NLOCATIONFILE+Last_Dash,".N.location");
*/
				break;
			default:
				usage();
				exit(0);
		}
	}	
	if (!Par_Count && !DICT) {usage();exit(1);}
	BP.Patternfile_Count --;
}

//}-----------------------------  Parse Command Line  -------------------------------------------------
void Read_INI(char* Config_File,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP)
{

	dictionary* Dictionary;
	char* Temp;
	char* Name;int Last_Dash,Temp_Int;char* Genome_Name;
	char* Command_Line_Buffer;

	if(!Config_File)
	{
		BP.Patternfile_Count= 0;
		BP.Misfile_Count= 0;
		BP.INSERTSIZE=200;
		BP.OFFSET=0;
		BP.PLUSSTRAND=FALSE;
		BP.LOADREVERSEONLY=FALSE;
		BP.SCANBOTH=FALSE;
		BP.ROLLOVER=FALSE;
		BP.MISHITS=FALSE;//TRUE;
		BP.ORIENTATIONS=2;//how many strand combinations to try..
		BP.USELOCATION = FALSE;
		BP.SMITH_WATERMAN = FALSE;
		BP.VERIFY=FALSE;//TRUE;//FALSE;
		BP.LOG=TRUE;
		BP.MAX_MISMATCHES= INT_MAX;
		BP.STD = 0;
		BP.SOLIDMAP=FALSE;
		BP.FLANKSIZE=0;
		BP.SW_FLANKSIZE=0;
		BP.FORCESOLID=0;
		BWAPOL=0;
		RECLIMIT=100000;
		UNIQ_HIT=FALSE;

		F.OUTPUTFILE=OUTFILE;
		Dictionary=iniparser_load("penguin.ini",FALSE);
		if (!Dictionary) Dictionary=iniparser_load("~/penguin.ini",FALSE);
	}
	else Dictionary=iniparser_load(Config_File,FALSE);

	if (Dictionary)
	{
		DICT=TRUE;
		if(Temp=iniparser_getstring(Dictionary,"files:head",NULL)) {BP.PATTERNFILE=Temp;BP.Patternfile_Count++;BP.PAIRING_MODE=NORMALFILEMODE;}
		if(Temp=iniparser_getstring(Dictionary,"files:tail",NULL)) {BP.PATTERNFILE1=Temp;BP.Patternfile_Count++;BP.PAIRING_MODE=NORMALFILEMODE;}
		if(Temp=iniparser_getstring(Dictionary,"files:genome",NULL)) 
		{
			Name=Temp;Last_Dash=0;Genome_Name=Temp;
			for(;Name[0]!=0;Name++)
			{
				if (Name[0]=='/') 
				{
					Last_Dash++;Genome_Name=Name;
				}
			}

			Command_Line_Buffer=(char*)malloc(6000);
			F.REVBWTINDEX = (char*)Command_Line_Buffer;
			if(Last_Dash) Last_Dash=Genome_Name-Temp+1; else Genome_Name--;
			strncpy(F.REVBWTINDEX,Temp,Last_Dash);
			F.REVBWTINDEX[Last_Dash+0]='r';F.REVBWTINDEX[Last_Dash+1]='e';F.REVBWTINDEX[Last_Dash+2]='v';
			strcpy(F.REVBWTINDEX+Last_Dash+3,Genome_Name+1);
			strcat(F.REVBWTINDEX+Last_Dash+3,".bwt"); 

			F.BWTFILE=F.REVBWTINDEX+500;
			strncpy(F.BWTFILE,Temp,Last_Dash);
			strcpy(F.BWTFILE+Last_Dash,Genome_Name+1);
			strcat(F.BWTFILE+Last_Dash,".bwt"); 


			F.REVOCCFILE = F.BWTFILE+500;
			strncpy(F.REVOCCFILE,Temp,Last_Dash);
			F.REVOCCFILE[Last_Dash+0]='r';F.REVOCCFILE[Last_Dash+1]='e';F.REVOCCFILE[Last_Dash+2]='v';
			strcpy(F.REVOCCFILE+Last_Dash+3,Genome_Name+1);
			strcat(F.REVOCCFILE+Last_Dash+3,".fmv"); 


			F.OCCFILE=F.REVOCCFILE+500;			
			strncpy(F.OCCFILE,Temp,Last_Dash);
			strcpy(F.OCCFILE+Last_Dash,Genome_Name+1);
			strcat(F.OCCFILE+Last_Dash,".fmv"); 

			F.SAFILE=F.OCCFILE+500;			
			strncpy(F.SAFILE,Temp,Last_Dash);
			strcpy(F.SAFILE+Last_Dash,Genome_Name+1);
			strcat(F.SAFILE+Last_Dash,".sa");

			F.REVSAFILE = F.SAFILE+500;
			strncpy(F.REVSAFILE,Temp,Last_Dash);
			F.REVSAFILE[Last_Dash+0]='r';F.REVSAFILE[Last_Dash+1]='e';F.REVSAFILE[Last_Dash+2]='v';
			strcpy(F.REVSAFILE+Last_Dash+3,Genome_Name+1);
			strcat(F.REVSAFILE+Last_Dash+3,".sa"); 

			F.BINFILE=F.REVSAFILE+500;			
			strncpy(F.BINFILE,Temp,Last_Dash);
			strcpy(F.BINFILE+Last_Dash,Genome_Name+1);
			strcat(F.BINFILE+Last_Dash,".pac");

			if(!BP.USELOCATION)
			{

				F.LOCATIONFILE=F.BINFILE+500;			
				strncpy(F.LOCATIONFILE,Temp,Last_Dash);
				strcpy(F.LOCATIONFILE+Last_Dash,Genome_Name+1);
				strcat(F.LOCATIONFILE+Last_Dash,".ann.location");
			}

			F.BLKFILE = F.LOCATIONFILE+500;
			strncpy(F.BLKFILE,Temp,Last_Dash);
			strcpy(F.BLKFILE+Last_Dash,Genome_Name+1);
			strcat(F.BLKFILE+Last_Dash,".blk.");

			F.INDFILE = F.BLKFILE+500;
			strncpy(F.INDFILE,Temp,Last_Dash);
			strcpy(F.INDFILE+Last_Dash,Genome_Name+1);
			strcat(F.INDFILE+Last_Dash,".ind.");

			F.RANGEFILE = F.INDFILE+500;
			strncpy(F.RANGEFILE,Temp,Last_Dash);
			strcpy(F.RANGEFILE+Last_Dash,Genome_Name+1);
			strcat(F.RANGEFILE+Last_Dash,".range");

			F.NLOCATIONFILE = F.RANGEFILE+500;
			strncpy(F.NLOCATIONFILE,Temp,Last_Dash);
			strcpy(F.NLOCATIONFILE+Last_Dash,Genome_Name+1);
			strcat(F.NLOCATIONFILE+Last_Dash,".N.location");
		}
		if(Temp=iniparser_getstring(Dictionary,"files:outputfile",NULL)) F.OUTPUTFILE=Temp; 
		if(Temp=iniparser_getstring(Dictionary,"files:logfile",NULL)) 
		{
			LOG_SUCCESS_FILE=Temp; 
			Log_SFile=File_Open(LOG_SUCCESS_FILE,"w");
		}
		if(Temp=iniparser_getstring(Dictionary,"files:location",NULL)) {F.LOCATIONFILE=Temp;BP.USELOCATION=TRUE;}
		if(Temp=iniparser_getstring(Dictionary,"files:headmis",NULL)) {BP.MISHITS=TRUE;BP.MISFILE1=Temp;BP.Misfile_Count++;} 
		if(Temp=iniparser_getstring(Dictionary,"files:tailmis",NULL)) {BP.MISHITS=TRUE;BP.MISFILE2=Temp;BP.Misfile_Count++;} 
		if(Temp=iniparser_getstring(Dictionary,"files:unmapped",NULL)) {WRITE_UNMAPPED=TRUE;BP.UNMAPPED_FILE=Temp;} 
		if(Temp=iniparser_getstring(Dictionary,"files:single",NULL)) {WRITE_SINGLE=TRUE;F.SINGLEFILE=Temp;} 
		if(Temp=iniparser_getstring(Dictionary,"files:multi",NULL)) {WRITE_MULTI=TRUE;BP.MULTI_FILE=Temp;} 
		if(Temp=iniparser_getstring(Dictionary,"files:disc",NULL)) {WRITE_DISCORDANT=TRUE;F.DISCORDANTFILE=Temp;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:insertsize",0)) {BP.INSERTSIZE=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:maxhits",0)) {MAXCOUNT=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:editdist",0)) {SW_EDIT_DIST_MAX=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:matchlen",0)) {SW_MIN_MATCH_LEN=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:ncount",0))	{NCOUNT=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:ntorand",0)){NPOLICY=true;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:nhits",0)) 
		{
			if((HITS_IN_SAM=Temp_Int)<=0) {printf ("Read_INI():Invalid parameter for -r\n");exit(-1);}
		}
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:verbose",0)) 
		{
			int Verb=atoi(optarg);
			if (Verb == 1){MISC_VERB=FALSE;} else if (Verb ==2){PROGRESSBAR=0;}
		}
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:maxtags",0)) {READS_TO_PROCESS=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:noheader",0)) {PRINT_DICTIONARY=FALSE;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:std",0)) {BP.STD=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:swon",0)) {BP.SMITH_WATERMAN=TRUE;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:flanksize",0)) {BP.SW_FLANKSIZE=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:forcelength",0)) {TRIM_LENGTH=Temp_Int;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:solidmap",0)) {BP.FORCESOLID=TRUE;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:bwapol",0)) {BWAPOL=TRUE;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:unique",0)) {UNIQ_HIT=TRUE;if (MAXCOUNT==1) MAXCOUNT=2;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:lowestmis",0)) {PRIORITYMODE=TRUE;} 
		if(Temp_Int=iniparser_getint(Dictionary,"parameters:maxmismatches",0)) 
		{
			BP.MAX_MISMATCHES =Temp_Int;
			if (BP.MAX_MISMATCHES <0 or BP.MAX_MISMATCHES >MAX_MISMATCHES_BOUND) BP.MAX_MISMATCHES=5;
		}
		
		if (BP.Patternfile_Count) BP.Patternfile_Count--;
	}
}
