#include <stdio.h>
#include <stdlib.h>
#include "bfix.h"
#include "rqindex.h"
#include <ctype.h>
#include "common.h"
extern int Genome_Position;
extern int Genome_Count;
extern unsigned Conversion_Factor;
extern char* LOG_SUCCESS_FILE;
extern FILE* Log_SFile;
unsigned long Block_Size = 0;
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Hits_File
 *  Description:  opens the BAT file containing hits and passed the information in it 
 *  		  to some global variables... else load mate files..
 * =====================================================================================
 */
void Load_Hits_File(BATPARAMETERS BP,FILELIST & FList,In_File & IN,char MODE)
{
#define READ_ASSERT(X,Y) {if ((X)<(Y)) {if(LOG_SUCCESS_FILE)fprintf(Log_SFile,"Load_Hits_File(): Read error...\n"); printf("Load_Hits_File(): Read error...\n");exit(-1);}}

if (BATFILEMODE == MODE)
{
	Header Head;
	char MAX_MISMATCHES,ROLLOVER,SCANBOTH;

	FList.Head=File_Open(BP.PATTERNFILE,"rb");
	IN.File_Size=Get_File_Size(FList.Head);
	READ_ASSERT(fread(&Head,1,sizeof(Header),FList.Head),sizeof(Header));
	if(!(Head.ID[0]=='B'&&Head.ID[1]=='A'&&Head.ID[2]=='T')) 
	{
		if (LOG_SUCCESS_FILE) fprintf(Log_SFile,"Not a BAT file\n");
		printf("Not a BAT file\n");
		exit(-1);
	}
	IN.MAXHITS=Head.MAXHITS;
	IN.STRINGLENGTH=Head.Tag_Length;
	IN.PRINT_DESC=Head.Print_Desc;
	IN.LOADREVERSEONLY=Head.Index_Count;

	//read(Data_File,&MAX_MISMATCHES,sizeof(MAX_MISMATCHES));if (MAX_MISMATCHES >5) Stat_Size=7*sizeof(unsigned short); else Stat_Size=(MAX_MISMATCHES+1)*sizeof(unsigned short);
	READ_ASSERT(fread(&MAX_MISMATCHES,sizeof(MAX_MISMATCHES),1,FList.Head),1);if (MAX_MISMATCHES >5) IN.Stat_Size=7*sizeof(unsigned short); else IN.Stat_Size=(MAX_MISMATCHES+1)*sizeof(unsigned short);
	//gzread(Data_File,&TAG_COPY_LEN,sizeof(int));
	READ_ASSERT(fread(&IN.TAG_COPY_LEN,sizeof(int),1,FList.Head),1);
	//gzread(Data_File,&ROLLOVER,sizeof(char));
	READ_ASSERT(fread(&ROLLOVER,sizeof(char),1,FList.Head),1);
	//gzread(Data_File,&SCANBOTH,sizeof(char));
	READ_ASSERT(fread(&SCANBOTH,sizeof(char),1,FList.Head),1);


	//gzread(Data_File,&Length_Array[1],sizeof(int));
	READ_ASSERT(fread(&IN.Length_Array[1],sizeof(int),1,FList.Head),1); IN.HEAD_LENGTH=IN.Length_Array[1];
	//gzread(Data_File,&Length_Array[2],sizeof(int));
	READ_ASSERT(fread(&IN.Length_Array[2],1,sizeof(int),FList.Head),sizeof(int)); IN.TAIL_LENGTH=IN.Length_Array[2];
	IN.FILETYPE=Head.FILETYPE;
}
else
{
	FList.Head=File_Open(BP.PATTERNFILE,"rb");
	FList.Tail=File_Open(BP.PATTERNFILE1,"rb");
}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Pair
 *  Description:  Read hits for Head and Tail of a tag into Head_HitsX and Tail_HitsX resply
 *  		  Terminates lists when first encountering Head/Tail_Hits[x].Start=0.
 *  		  returns true if hits read...
 * =====================================================================================
 */

#define REPEAT 1
#define PM 2
#define PP 3
#define MP 4
#define MM 5

char Read_Pair(BWT* fwfmi,BWT* revfmi,SARange* Head_Hits_Pos,SARange* Head_Hits_Neg, SARange* Tail_Hits_Pos,SARange* Tail_Hits_Neg,FILE* Data_File,In_File & IN,Record_Info & Hit_Details,char & MAX_HIT_FAULT,unsigned Conversion_Factor)
{
#undef READ_ASSERT
#define READ_ASSERT(X,Y) {if ((X)<(Y)) {if(LOG_SUCCESS_FILE)fprintf(Log_SFile,"Read_Pair(): Read error...\n"); printf("Read_Pair(): Read error...\n");exit(-1);}}
	MAX_HIT_FAULT=FALSE;
	int Head_Hit_Pos_Count=0;
	int Head_Hit_Neg_Count=0;
	int Tail_Hit_Pos_Count=0;
	int Tail_Hit_Neg_Count=0;
	unsigned i, Bitsread;
	unsigned Start,End;
	unsigned Location;
	int Desc_End;
	unsigned Previous_Tag,Hits;
	char letter;
	unsigned pos;
	int StringLength;
	unsigned Progress=0,Average_Length;
	unsigned Number_of_Tags=100;
	char Mismatches_Desc[1000];
	char Quality[MAXTAG+1];
	char* Mismatch_Desc_Ptr;
	char Last_Orientation,Last_Half;
	char* Ins_Format;
	char* Del_Format;
	char* Mis_Format;
	static char First_Pass=TRUE;
	static int Break=0;
	char Tag_Type;
	char* Seek_Tail;
	char N[500];
	unsigned short Stats1[7];
	unsigned short Stats2[7];
	char New_Record[4];
	Mismatches_Record_GIS MismatchesGIS;
	Output_Record Record;

	if (First_Pass) {READ_ASSERT(fread(&Tag_Type,1,1,Data_File),1);} 
	if(fgets(Hit_Details.Description,MAXDES,Data_File)==NULL) {if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Read_Pair(): error reading file...\n");printf("Read_Pair(): error reading file...\n");exit(-1);};
	Seek_Tail=Hit_Details.Description;
	while(*(Seek_Tail)!='\n'&&*(Seek_Tail)!='\t'&&*(Seek_Tail)!=' ') Seek_Tail++;*Seek_Tail=0;//make description to a string...
	READ_ASSERT(fread(Stats1,IN.Stat_Size,1,Data_File),1);
	READ_ASSERT(fread(Stats2,IN.Stat_Size,1,Data_File),1);
	READ_ASSERT(fread(IN.Tag_Copy,IN.TAG_COPY_LEN,1,Data_File),1);//if(NORMAL_TAGS && FILETYPE==FQ) gzread(Data_File,Quality,TAG_COPY_LEN);
	if(IN.FILETYPE==FQ) READ_ASSERT(fread(Hit_Details.Quality,1,IN.TAG_COPY_LEN,Data_File),IN.TAG_COPY_LEN);
	READ_ASSERT(fread(&Tag_Type,1,1,Data_File),1); 
	if (First_Pass)
	{
		IN.Positive_Tail=IN.Tag_Copy;First_Pass=FALSE;
		while (*(IN.Positive_Tail++)!='\t'){Break++;};
	}
	IN.Tag_Copy[Break]=0;


//--------------------------------------------------------------------------------------------
	for(;;)
	{
		//fread(&Tag_Type,1,1,Data_File);
		if('&' == Tag_Type) return FALSE;
		if('@' == Tag_Type)//start of new tag
		{
			Head_Hits_Pos[Head_Hit_Pos_Count].Start=0;
			Head_Hits_Neg[Head_Hit_Neg_Count].Start=0;
			Tail_Hits_Pos[Tail_Hit_Pos_Count].Start=0;
			Tail_Hits_Neg[Tail_Hit_Neg_Count].Start=0;
			if((Head_Hit_Pos_Count + Head_Hit_Neg_Count ==1) && (Tail_Hit_Pos_Count + Tail_Hit_Neg_Count==1)) //unique pair..
			{
				if (Head_Hit_Pos_Count ==1 && !Head_Hit_Neg_Count)// H + uniq
				{
					if (Tail_Hit_Pos_Count ==1 && !Tail_Hit_Neg_Count)//T + uniq
					{
						return PP;
					}
					else if (Tail_Hit_Pos_Count ==0 && Tail_Hit_Neg_Count==1)//T - uniq
					{
						return PM;
					}
				}
				else if (Head_Hit_Pos_Count ==0 && Head_Hit_Neg_Count ==1)// H - uniq
				{
					if (Tail_Hit_Pos_Count ==1 && !Tail_Hit_Neg_Count)//T + uniq
					{
						return MP;
					}
					else if (Tail_Hit_Pos_Count ==0 && Tail_Hit_Neg_Count==1)//T - uniq
					{
						return MM;
					}
				}
			}
			return REPEAT;
		}		
		else//Process hit....
		{
			if ('&' == Tag_Type) return FALSE;
//--------------------------------------------------------------------------------------------
			//gzread(Data_File,New_Record,3);
			READ_ASSERT(fread(New_Record,1,3,Data_File),3);
			//if (New_Record[2]) gzread(Data_File,N,New_Record[2]*2);
			if (New_Record[2]) READ_ASSERT(fread(N,1,New_Record[2]*2,Data_File),New_Record[2]*2);
			//gzread(Data_File,&Record,sizeof(Output_Record));
			READ_ASSERT(fread(&Record,1,sizeof(Output_Record),Data_File),sizeof(Output_Record));
			//gzread(Data_File,&MismatchesGIS,Record.Mismatches+sizeof(unsigned));
			READ_ASSERT(fread(&MismatchesGIS,1,Record.Mismatches+sizeof(unsigned),Data_File),Record.Mismatches+sizeof(unsigned));
			StringLength=IN.Length_Array[New_Record[0]];
			
			if (!Record.Gap)
			{
				if (Record.Skip) Record.Start = Conversion_Factor-revfmi->saValue[Record.Start/revfmi->saInterval]+Record.Skip-1;
				else Record.Start=Conversion_Factor-BWTSaValue(revfmi,Record.Start);
			}
			if(New_Record[0]==1)
			{
				if (New_Record[1] == '-')
				{
					Head_Hits_Neg[Head_Hit_Neg_Count].Start=Record.Start;
					Head_Hits_Neg[Head_Hit_Neg_Count].End=Record.Start+Record.Gap;
					Head_Hits_Neg[Head_Hit_Neg_Count].Mismatches=Record.Mismatches;
					Head_Hit_Neg_Count++;
				}
				else
				{
					Head_Hits_Pos[Head_Hit_Pos_Count].Start=Record.Start;
					Head_Hits_Pos[Head_Hit_Pos_Count].End=Record.Start+Record.Gap;
					Head_Hits_Pos[Head_Hit_Pos_Count].Mismatches=Record.Mismatches;
					Head_Hit_Pos_Count++;
				}
#ifdef DEBUG
				if(Head_Hit_Neg_Count >= MAX_HITS_TO_STORE) {Head_Hit_Neg_Count--;MAX_HIT_FAULT=TRUE;}else if( Head_Hit_Pos_Count >= MAX_HITS_TO_STORE) {Head_Hit_Pos_Count--;MAX_HIT_FAULT=TRUE;}//{printf("Read_Pair: Too many hits !..\n");exit(1);}
#endif
			}
			else
			{
				if (New_Record[1] == '-')
				{
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Start=Record.Start;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].End=Record.Start+Record.Gap;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Mismatches=Record.Mismatches;
					Tail_Hit_Neg_Count++;
				}
				else
				{
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Start=Record.Start;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].End=Record.Start+Record.Gap;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Mismatches=Record.Mismatches;
					Tail_Hit_Pos_Count++;
				}
#ifdef DEBUG
				if(Tail_Hit_Neg_Count >= MAX_HITS_TO_STORE){Tail_Hit_Neg_Count--;MAX_HIT_FAULT=TRUE;} else if( Tail_Hit_Pos_Count >= MAX_HITS_TO_STORE) {Tail_Hit_Pos_Count--;MAX_HIT_FAULT=TRUE;}//printf("Read_Pair: Too many hits !..\n");exit(1);}
#endif
			}
		}
		READ_ASSERT(fread(&Tag_Type,1,1,Data_File),1);

	}
}


void Load_Range_Index(RQINDEX & R,int STRINGLENGTH,FMFILES F,unsigned & Entries)
{
#undef READ_ASSERT
#define READ_ASSERT(X,Y) {if ((X)<(Y)) {if(LOG_SUCCESS_FILE)fprintf(Log_SFile,"Load_Range_Index(): Read error...\n"); printf("Load_Range_Index(): Read error...\n");exit(-1);}}
	unsigned long Index_Size = 0; //,Block_Size;
	int T1=strlen(F.INDFILE);
	int T2=strlen(F.BLKFILE);

	sprintf(F.INDFILE+strlen(F.INDFILE),"%d",STRINGLENGTH);
	sprintf(F.BLKFILE+strlen(F.BLKFILE),"%d",STRINGLENGTH);
	FILE* Index=File_Open(F.INDFILE,"rb");
	FILE* Blocks=File_Open(F.BLKFILE,"rb");
	R.SA_Index=(SA*)malloc((Index_Size=Get_File_Size(Index)));//contains the index to blocks...
	R.SA_Blocks=(int*)malloc((Block_Size=Get_File_Size(Blocks)));
	//memset(R.SA_Blocks, 0, Block_Size);
	if (!R.SA_Index || !R.SA_Blocks)
	{
		if(LOG_SUCCESS_FILE) fprintf(Log_SFile,"Load_Range_Index(): Memory allocation failed!..\n");
		printf ("Load_Range_Index(): Memory allocation failed!..\n");
		exit(-1);
	}
	READ_ASSERT(fread(&R.COMPRESS,1,1,Blocks),1);
	READ_ASSERT(fread(&Entries,sizeof(Entries),1,Blocks),1);
	READ_ASSERT(fread(R.SA_Index,Index_Size,1,Index),1);
	READ_ASSERT(fread(R.SA_Blocks,Block_Size-sizeof(Entries)-1,1,Blocks),1);

	*(F.INDFILE+T1)=0;
	*(F.BLKFILE+T2)=0;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Small_Gap
 *  Description:  Do the searching when at least one of the pairs have small SA range... 
 * 		  Store the result in Pairs, starting from Pairs[Pairs_Index]
 * 		  Terminates hit list when encountering Pairs[x].Head=0... 
 * =====================================================================================
 */
void Search_Small_Gap(BWT* revfmi,RQINDEX & R, SARange & Head,  SARange & Tail, int d,Pair* Pairs,int & Pairs_Index,unsigned MAXCOUNT,int & HITS,unsigned Entries,unsigned Conversion_Factor)
{
	unsigned H1,H2,T1,T2,L,H,M;
	Tag_Info Head_Info,Tail_Info;


	if(Head.Start==Head.End)//Head is unique...
	{

		H1=Head.Start;//Conversion_Factor-BWTSaValue(revfmi,Head.Start);
		if (Tail.Start==Tail.End)// both hits are unique...
		{
			T2=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
			if (T2>H1 && T2 < H1+d)//modify for multiple hits...
			{
				Pairs[Pairs_Index].Head=H1;
				Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
				Pairs[Pairs_Index].Tail=T2;
				Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
				Pairs_Index++;HITS++;
#ifdef DEBUG
				if(H1>T2) 
				{
					printf("Search_Small_Gap(6):Enum error...\n");
					exit(0);
				}
#endif
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else //tail has multiples...
		{
			if (Tail.End-Tail.Start > SAGAP_CUTOFF)
			{
				Load_Info(R,Tail_Info,Tail,Entries);
				T1=Tail_Info.First;T2=Tail_Info.Last;


				if(H1<T1 )//Possible case for T1>H1 
				{
					if(H1+d>T1)
					{
						M=0;
						while (H1<T1 && H1+d>=T1)//enumerate hits...
						{
							Pairs[Pairs_Index].Head=H1;
							Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
							Pairs[Pairs_Index].Tail=T1;
							Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
							Pairs_Index++;HITS++;
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(5):Enum error...\n");
							exit(0);
						}
#endif
							if (HITS>=MAXCOUNT) 
							{
								Pairs[Pairs_Index].Head=0;
								return;
							}
							M++;
							T1=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);
							if(M>=Tail_Info.Gap) break;
						}
					}
				}
				else if(T2>H1) //H1 inside tail gaps..
				{

					L=0;
					H=Tail_Info.Gap;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor) > H1)
						{
							H=M;
						}
						else
						{
							L=M+1;
						}
					}
					if (L==H) M=H;//find tail position closest to unique head...

					T1=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);
					while (H1<T1 && H1+d>=T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
						Pairs_Index++;HITS++;
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(4):Enum error...\n");
							exit(0);
						}
#endif
						if (HITS>=MAXCOUNT) 
						{
							Pairs[Pairs_Index].Head=0;
							return;
						}
						M++;
						T1=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);
						if(M>=Tail_Info.Gap) break;
					}
				}

				Pairs[Pairs_Index].Head=0;
				return;
			}
			else//Unique head and tail with gap below cutoff...
			{
				for (unsigned i=Tail.Start;i<=Tail.End;i++)
				{
					unsigned Hit=(Conversion_Factor-BWTSaValue(revfmi,i));
					if(Hit > H1 && Hit < H1+d)
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
						Pairs[Pairs_Index].Tail=Hit;
						Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
						Pairs_Index++;HITS++;
#ifdef DEBUG
						if(H1>Hit) 
						{
							printf("Search_Small_Gap(3):Enum error...\n");
							exit(0);
						}
#endif

						if(HITS >=MAXCOUNT) break;
					}
				}
				Pairs[Pairs_Index].Head=0;
				return;
			}
		}
	}
	else //if(Tail.End==Tail.Start)//Tail is unique...
	{
		T1=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
		if(Head.End-Head.Start>SAGAP_CUTOFF)//Unique tail, but with multiple possible heads...
		{
			Load_Info(R,Head_Info,Head,Entries);
			H1=Head_Info.First;H2=Head_Info.Last;
			if(T1 > H1) //Head should not be after T1...
			{
				if(H2+d >T1)//Tail not too far away from heads...
				{
					//T1-d is between H1 and H2, search for the closest hit...
					L=0;
					H=Head_Info.Gap-1;
					unsigned T1_Temp=T1;
					if (T1>d) T1-=d; else T1=0;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(revfmi,R,Head_Info,M,Conversion_Factor) < T1)
						{
							L=M+1;
						}
						else
						{
							H=M;
						}
					}
					if (L==H) M=L;
					H1=Get_Location(revfmi,R,Head_Info,M,Conversion_Factor);

					//T1+=d;
					T1=T1_Temp;
					while (H1<T1 && H1+d>=T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
						Pairs_Index++;HITS++;
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(2):Enum error...\n");
							exit(0);
						}
#endif
						if (HITS>=MAXCOUNT) 
						{
							break;	
						}
						M++;
						H1=Get_Location(revfmi,R,Head_Info,M,Conversion_Factor);
						if(M>=Head_Info.Gap) break;
					}

				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else//Unique tail and head with gap below cutoff...
		{
			for (unsigned i=Head.Start;i<=Head.End;i++)//try all hits...
			{
				unsigned Hit=(Conversion_Factor-BWTSaValue(revfmi,i));
				if(Hit < T1 && T1 <= Hit+d)
				{
					Pairs[Pairs_Index].Head=Hit;
					Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
					Pairs[Pairs_Index].Tail=T1;
					Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
					Pairs_Index++;HITS++;
#ifdef DEBUG
				if(Hit>T1) 
				{
					printf("Search_Small_Gap(1):Enum error...\n");
					exit(0);
				}
#endif
					if(HITS >=MAXCOUNT) break;
				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Head_Tail...
 *  Description:  find the possible pairs for head and tail that are a distance d apart.. 
 *  		  Pairs are stored in the structure Pair..
 *  		  Pairs_Index is the offset of the sentinel indicating last hit...
 * =====================================================================================
 */


void Get_Head_Tail(BWT* revfmi,RQINDEX & R,SARange & Head, SARange & Tail,int d,Pair* Pairs,int & Pairs_Index,unsigned MAXCOUNT,int & HITS,unsigned Entries,unsigned Conversion_Factor)
{
	unsigned H1,H2,T1,T2;
	unsigned H,M,L;
	unsigned Location=0,Hd,Tl;

	Tag_Info Head_Info,Tail_Info;


	Load_Info(R,Head_Info,Head,Entries);
	Load_Info(R,Tail_Info,Tail,Entries);
	H1=Head_Info.First;H2=Head_Info.Last + d;
	T1=Tail_Info.First;T2=Tail_Info.Last;


	if (T1>H2 || T2<H1) //disjoint head and pair...
	{
		Pairs[Pairs_Index].Head=0;//cannot be paired...
		return;
	}

	if (T1>H1) M=0;//Get in M the index of the closest value of Tail larger than H1
	else
	{
		L=0;
		H=Tail_Info.Gap-1;
		while (L < H)
		{
			M=(L+H)/2;
			if (Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor) < H1)
			{
				L=M+1;
			}
			else
			{
				H=M;
			}
		}
		if (L==H) M=H;
	}


	Hd=H1;H=1;
	Tl=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);M++;

#ifdef DEBUG
	if(Tl<Hd) {printf("Get_Head_Tail(): bin search error..\n");exit(0);}
#endif

	int MHead=1;//Array pos. of head...
	int TempM;
	unsigned TempTl;

	while (Hd<T2)//Tail must always be after Head...
	{
		if (Hd<Tl)
		{
			TempM=M;TempTl=Tl;
			while (Hd+d>=Tl)//enum. Hits...
			{
				Pairs[Pairs_Index].Head=Hd;
				Pairs[Pairs_Index].MismatchesH=Head.Mismatches;
				Pairs[Pairs_Index].Tail=Tl;
				Pairs[Pairs_Index].MismatchesT=Tail.Mismatches;
				Pairs_Index++;HITS++;
#ifdef DEBUG
				if(Hd>Tl) 
				{
					return;
					//printf("Get_Head_Tail():Enum error...\n");
					exit(0);
				}
#endif
				if (HITS>=MAXCOUNT) 
				{
					Pairs[Pairs_Index].Head=0;
					return;
				}
				if(M >= Tail_Info.Gap) break;else Tl=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);
				M++;
			}
			M=TempM;Tl=TempTl;

			if (MHead>=Head_Info.Gap) break; else Hd=Get_Location(revfmi,R,Head_Info,MHead,Conversion_Factor);
			MHead++;
		}
		else
		{
			if(M >= Tail_Info.Gap) break;else Tl=Get_Location(revfmi,R,Tail_Info,M,Conversion_Factor);
			M++;
		}
	}
	Pairs[Pairs_Index].Head=0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Info
 *  Description:  Loads information necessary for reading coded info.... 
 * =====================================================================================
 */

void Load_Info( RQINDEX & R,Tag_Info & Tag, SARange & Head,unsigned Entries)
{
	Tag.Gap=Head.End-Head.Start+1;
	Tag.SA_Start=Head.Start;
	Tag.Block_Start=Get_Block_Start(R,Head.Start, Tag.Index,Entries);//get info block..
	Tag.First=R.SA_Index[Tag.Index].Start_Location;
	Tag.Last=R.SA_Index[Tag.Index].End_Location;
	Tag.Field_Length=log2(Tag.Gap);
}


unsigned Get_Location(BWT* revfmi,RQINDEX & R,Tag_Info & Tag, unsigned Offset,unsigned Conversion_Factor)
{
	unsigned SAPos;
	if(0==Offset) return Tag.First;
	if (Offset==Tag.Gap-1) return Tag.Last;
	Offset--;
	if(R.COMPRESS)
	{
		SAPos=Tag.SA_Start + bfx((unsigned char*)R.SA_Blocks,Tag.Block_Start+(Offset*Tag.Field_Length),Tag.Field_Length);
		return Conversion_Factor-BWTSaValue(revfmi,SAPos);
	}
	else
	{
		//if(Tag.Block_Start+Offset >= Block_Size/sizeof(int) && Tag.Block_Start+Offset < Block_Size) printf("\nSIZE %ld %ld %ld %ld %d %p\n", Tag.Block_Start, Offset, Block_Size, Tag.Block_Start+Offset, Tag.First, R.SA_Blocks);
		if(Tag.Block_Start+Offset >= Block_Size) return Tag.First;
		return (unsigned)R.SA_Blocks[Tag.Block_Start+Offset];
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Block_Start
 *  Description:  Finds the Block Value of the SAValue by bin search...
 *  		  Assumes SAValue exists...
 * =====================================================================================
 */
unsigned Get_Block_Start(RQINDEX & R,unsigned SAValue,unsigned & M,unsigned H)
{
	unsigned L=0;

	while (L < H)
	{
		M=(L+H)/2;
		if (R.SA_Index[M].Start < SAValue)
		{
			L=M+1;
		}
		else
		{
			H=M;
		}
	}

	if (L==H) M=H; 

#ifdef DEBUG
	if (R.SA_Index[M].Start!=SAValue) 
	{
		return R.SA_Index[M].End;
		//printf("Get_Block_Start(): Bad SAValue !\n");
		exit(0);
	}
#endif

	return R.SA_Index[M].End;
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_PairV
 *  Description:  Read hits for Head and Tail of a tag into Head_Hits and Tail_Hits resply
 *  		  Terminates lists when first encountering Head/Tail_Hits[x].Start=0.
 *  		  returns true if hits read...
 * =====================================================================================
 */
char Read_PairV(SARange* Head_Hits_Pos,SARange* Head_Hits_Neg, SARange* Tail_Hits_Pos,SARange* Tail_Hits_Neg,FILE* Data_File,In_File IN,Record_Info & Hit_Details)
{
	char Layout[2];//head/tail, orientation...
	char Tag_Type;
	int Head_Hit_Neg_Count=0;
	int Head_Hit_Pos_Count=0;
	int Tail_Hit_Neg_Count=0;
	int Tail_Hit_Pos_Count=0;
	int Desc_End;
	char Mismatch_Count;
	//char Description[MAXDES+1];
	char* Dest_String;
	char* Source_String;
	static char First_Pass=TRUE;
	Mismatches_Record Mismatches;
	Output_Record Record;
	char Original_String[MAXSTRINGLENGTH+1];
	
	for(;;)	
	{

		READ_ASSERT(fread(&Tag_Type,1,1,Data_File),1);
		if('&' == Tag_Type) break;
		if('@' == Tag_Type)//start of new tag
		{
			
			//Pass=1;PM='+';
			if (First_Pass) 
			{
				if(IN.PRINT_DESC) 
				{
					if(!fgets(Hit_Details.Description,MAXDES+1,Data_File)) {printf("Read_PairV():Error reading file..\n");exit(0);};
					for(Desc_End=0;Hit_Details.Description[Desc_End]!=' ' && Hit_Details.Description[Desc_End]!='\t' && Hit_Details.Description[Desc_End]!='\n' && Hit_Details.Description[Desc_End]!='\r' && Hit_Details.Description[Desc_End]!=0;Desc_End++);
					Hit_Details.Description[Desc_End]=0;
				};

				if(!fgets(Original_String,MAXSTRINGLENGTH,Data_File)){printf("Read_PairV(): Error reading file..\n");exit(0);};
				for(char* Next_String=Original_String;Next_String[0]!='\n' && Next_String[0]!='\r' && Next_String[0]!=0;Next_String++) Next_String[0]=tolower(Next_String[0]);
				memcpy(Hit_Details.Positive_Head,Original_String,IN.HEAD_LENGTH);
				memcpy(Hit_Details.Positive_Tail,Original_String+IN.HEAD_LENGTH+1,IN.TAIL_LENGTH);
				/*for (int i=0;i<=HEAD_LENGTH-1;i++)
				{
					Hit_Details.Reverse_Head[HEAD_LENGTH-1-i]=Code_To_Char[Char_To_CodeC[Hit_Details.Positive_Head[i]]];
				}
				for (int i=0;i<=TAIL_LENGTH-1;i++)
				{
					Hit_Details.Reverse_Tail[TAIL_LENGTH-1-i]=Code_To_Char[Char_To_CodeC[Hit_Details.Positive_Tail[i]]];
				}*/

				First_Pass=FALSE;
				continue;
			} 
			else 
			{
				Head_Hits_Pos[Head_Hit_Pos_Count].Start=0;
				Head_Hits_Neg[Head_Hit_Neg_Count].Start=0;
				Tail_Hits_Pos[Tail_Hit_Pos_Count].Start=0;
				Tail_Hits_Neg[Tail_Hit_Neg_Count].Start=0;
				fseek(Data_File,-1,SEEK_CUR);
				First_Pass=TRUE;
				return TRUE;
			}
		}		
		else 
		{
			if ('&' == Tag_Type) return FALSE;
			READ_ASSERT(fread(&Layout,2,1,Data_File),1);
			/*if(Layout[0]=='H')//Head..
			{
				Dest_String=Head_String;STRINGLENGTH=HEAD_LENGTH;
				Source_String=Original_String;
			} 
			else//Tail..
			{
				Dest_String=Tail_String;STRINGLENGTH=TAIL_LENGTH;
				Source_String=Original_String+TAIL_LENGTH+1;
			}
			if (Layout[1] == '-')
			{
				for (unsigned i=0;i<=STRINGLENGTH-1;i++)
				{
					Dest_String[STRINGLENGTH-1-i]=Code_To_Char[Char_To_CodeC[Source_String[i]]];
				}
			}
			else
			{
				memcpy(Dest_String,Source_String,MAXSTRINGLENGTH);
			}*/

			READ_ASSERT(fread(&Mismatches,sizeof(Mismatches_Record),1,Data_File),1);
			READ_ASSERT(fread(&Record,sizeof(Output_Record),1,Data_File),1);
			Mismatch_Count = Record.Mismatches;

			/*if (Mismatch_Count)//mismatch positions..
			{
				if(Layout[1]=='+')
				{
					for( int i=0;i<Mismatch_Count;i++)
					{
						Dest_String[63 & (Mismatches.Mismatch_Pos>>(6*i))]=Code_To_Char[Mismatches.Mismatch_Char>>(2*i) & 3];
					}
				}
				else
				{
					for( int i=0;i<Mismatch_Count;i++)
					{
						Dest_String[(63 & (Mismatches.Mismatch_Pos>>(6*i)))]=Code_To_Char[Mismatches.Mismatch_Char>>(2*i) & 3];
					}
				}
			}*/

			if(Layout[0]=='H')
			{
				if (Layout[1] == '-')
				{
					Head_Hits_Neg[Head_Hit_Neg_Count].Start=Record.Start;
					Head_Hits_Neg[Head_Hit_Neg_Count].End=Record.Start+Mismatches.Gap;
					Head_Hits_Neg[Head_Hit_Neg_Count].Mismatches=Mismatch_Count;
					Head_Hits_Neg[Head_Hit_Neg_Count].Mismatch_PosX=Mismatches.Mismatch_PosX;
					Head_Hits_Neg[Head_Hit_Neg_Count].Mismatch_Char=Mismatches.Mismatch_Char;
					Head_Hit_Neg_Count++;
				}
				else
				{
					Head_Hits_Pos[Head_Hit_Pos_Count].Start=Record.Start;
					Head_Hits_Pos[Head_Hit_Pos_Count].End=Record.Start+Mismatches.Gap;
					Head_Hits_Pos[Head_Hit_Pos_Count].Mismatches=Mismatch_Count;
					Head_Hits_Pos[Head_Hit_Pos_Count].Mismatch_PosX=Mismatches.Mismatch_PosX;
					Head_Hits_Pos[Head_Hit_Pos_Count].Mismatch_Char=Mismatches.Mismatch_Char;
					Head_Hit_Pos_Count++;
				}
#ifdef DEBUG
				if(Head_Hit_Neg_Count >= MAX_HITS_TO_STORE || Head_Hit_Pos_Count >= MAX_HITS_TO_STORE) {printf("Read_PairV: Too many hits !..\n");exit(1);}
#endif
			}
			else
			{
				if (Layout[1] == '-')
				{
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Start=Record.Start;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].End=Record.Start+Mismatches.Gap;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Mismatches=Mismatch_Count;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Mismatch_PosX=Mismatches.Mismatch_PosX;
					Tail_Hits_Neg[Tail_Hit_Neg_Count].Mismatch_Char=Mismatches.Mismatch_Char;
					Tail_Hit_Neg_Count++;
				}
				else
				{
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Start=Record.Start;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].End=Record.Start+Mismatches.Gap;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Mismatches=Mismatch_Count;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Mismatch_PosX=Mismatches.Mismatch_PosX;
					Tail_Hits_Pos[Tail_Hit_Pos_Count].Mismatch_Char=Mismatches.Mismatch_Char;
					Tail_Hit_Pos_Count++;
				}
#ifdef DEBUG
				if(Tail_Hit_Neg_Count >= MAX_HITS_TO_STORE || Tail_Hit_Pos_Count >= MAX_HITS_TO_STORE) {printf("Read_PairV: Too many hits !..\n");exit(1);}
#endif
			}


			
			/*for (int j=0;j<=Mismatches.Gap && Hits< MAXHITS;j++)
			{
				if( Record.Index)//print record...
				{
					Location=Conversion_Factor-BWTSaValue(revfmi,Record.Start)-Offset;
				}
				else
				{
					Location=BWTSaValue(fmi,Record.Start)-Offset;
				}
				fseeko64(Binary_File,Location,SEEK_SET);
				fread(&Test_String,StringLength,1,Binary_File);
				if (strncmp(Test_String,String,StringLength))
				{
					for (unsigned i=0;i<=StringLength-1;i++)
					{
						String[i]=Code_To_Char[Char_To_CodeC[String[i]]];
					}

					if (strncmp(Test_String,RevString,StringLength))
					{
						printf("Error in Tag %d\t%s \n: Location %u is\t%s\t%s\n",Record.Tag, String,Location,Test_String,Description);
					}
				}
				Hits++;Record.Start++;Total_Hits++;
			}*/


		}
	}
printf("eek");
	return FALSE;
}
