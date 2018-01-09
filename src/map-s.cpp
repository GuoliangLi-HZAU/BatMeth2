#include <map>
#include "batlib-s.h"
/*
void Zig_Zag()
{
		if(Zero_Mismatch(B.Forward,L,revfmi,MF))
		{
			Mapped++;
		}
		else if(Zero_Mismatch(B.Complement,L,revfmi,MC))
		{
			Mapped++;
		}
		else if (One_Mismatch(B.Forward,L,MAXHITS,fwfmi,revfmi,MF))
		{
			Mapped++;
		}
		else if (One_Mismatch(B.Complement,L,MAXHITS,fwfmi,revfmi,MC))
		{
			Mapped++;
		}
		else if (Two_Mismatch(B.Forward,L,MAXHITS,fwfmi,revfmi,MF))
		{
			Mapped++;
		}
		else if (Two_Mismatch(B.Complement,L,MAXHITS,fwfmi,revfmi,MC))
		{
			Mapped++;
		}
		else if (Three_Mismatch(B.Forward,L,MAXHITS,fwfmi,revfmi,MF))
		{
			Mapped++;
		}
		else if (Three_Mismatch(B.Complement,L,MAXHITS,fwfmi,revfmi,MC))
		{
			Mapped++;
		}
		else if (Four_Mismatch(B.Forward,L,MAXHITS,fwfmi,revfmi,MF))
		{
			Mapped++;
		}
		else if (Four_Mismatch(B.Complement,L,MAXHITS,fwfmi,revfmi,MC))
		{
			Mapped++;
		}
		else 
		{
			printf ("%s\n",R.Tag_Copy);
		}
}
*/
/*void *Map_t(void* TP)
{
	unsigned Actual_Tag=0;
	unsigned Mapped=0;
	int LOOKUPSIZE;
	READ R;BATREAD B;
	GUESS G;
	MEMX MF,MC;
	Thread_Arg* T=(Thread_Arg*)TP;

//----------- copy arguments ------------------
	BWT* fwfmi=T->fwfmi;
	BWT* revfmi=T->revfmi;
	//FILE* Input_File=T->Input_File;
	FILE* Input_File=(T->In).Input_File;
	FILE* Output_File=T->Output_File;
	LEN L=T->L; 
	INFILE F=T->In;
	unsigned MAXHITS=(T->Bat_Para).MAXHITS;
	char ONEFMINDEX=(T->Bat_Para).ONEFMINDEX;
	char MAX_MISMATCHES=(T->Bat_Para).MAX_MISMATCHES;
	char SOLID=(T->In).SOLID; 
	MEMLOOK MLook=T->MLook;
	//int STRINGLENGTH=T->STRINGLENGTH;
	OUTPUT O=T->Output;

//----------- calculate read portions ------------------
	LOOKUPSIZE=Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	MF.Lookupsize=LOOKUPSIZE;MC.Lookupsize=LOOKUPSIZE;
	F.Buffer=(char*) malloc(6000);
	char* Output;unsigned Progress=0;unsigned Number_of_Tags=1000;
//----------- initialise memory structures ------------------
	Copy_MEM(MLook,MF,MC,MAX_MISMATCHES);
	if (T->ThreadID ==0) Init_Progress();

	while (Read_Tag(Input_File,FQ,R))
	{
		//if (ftello64(Input_File) >= Last) break;


		R.Tag_Number=1;R.Read_Number=Actual_Tag;
		Process_Read(R,B,MF,MC);
		Actual_Tag++;
		int Hits=0;
		Output=O.Buffer;
		if (T->ThreadID==0)
		{
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(F.Input_File);
				off64_t Average_Length=Current_Pos/Actual_Tag+1;//+1 avoids divide by zero..
				Number_of_Tags=(F.File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/F.File_Size);
			}
		}

		if(Hits=Guess_Orientation(fwfmi,revfmi,MF,MC,L,G,B))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
			if (MAXHITS <=Hits || MAX_MISMATCHES == 0) continue;
		}


		if(Hits=Map_Strand_Guess(MAX_MISMATCHES,MAXHITS,L,fwfmi,revfmi,G.Guessed))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
		}
		else if(Hits=Map_Strand_Guess(MAX_MISMATCHES,MAXHITS,L,fwfmi,revfmi,G.Guess_Complement))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guess_Complement,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guess_Complement,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
		}
		else
		{
			//printf ("%u\t%s\n",Actual_Tag,R.Tag_Copy);
		}


	}
	char End_Mark='&';
	fwriteX(&End_Mark,1,1,Output_File);//write sentinel..
	return (void*)(Mapped);
*/

/*unsigned Map(BWT* fwfmi,BWT* revfmi,FILE* Output_File,LEN & L, BATPARAMETERS P, MEMLOOK MLook,INFILE & F)
{
	unsigned Actual_Tag=0;
	unsigned Mapped=0;
	int LOOKUPSIZE;
	READ R;BATREAD B;
	GUESS G;
	MEMX MF,MC;

	unsigned MAXHITS=P.MAXHITS;
	char ONEFMINDEX=P.ONEFMINDEX;
	int MAX_MISMATCHES=P.MAX_MISMATCHES;

	F.Buffer=(char*) malloc(6000);
	char* Output;

//----------- calculate read portions ------------------
	LOOKUPSIZE=Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	MF.Lookupsize=LOOKUPSIZE;MC.Lookupsize=LOOKUPSIZE;
//----------- initialise memory structures ------------------
	Copy_MEM(MLook,MF,MC,MAX_MISMATCHES);
	int Progress=0;unsigned Number_of_Tags=1000;
	Init_Progress();

	while (Read_Tag(F.Input_File,F.FILETYPE,R))
	{
		//if (ftello64(F.Input_File) >= Last) break;
		Progress++;
		if (Progress>Number_of_Tags) 
		{
			off64_t Current_Pos=ftello64(F.Input_File);
			off64_t Average_Length=Current_Pos/Actual_Tag+1;//+1 avoids divide by zero..
			Number_of_Tags=(F.File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/F.File_Size);
		}
		R.Tag_Number=1;R.Read_Number=Actual_Tag;
		Process_Read(R,B,MF,MC);
		Actual_Tag++;
		int Hits=0;Output=F.Buffer;

		if(Hits=Guess_Orientation(fwfmi,revfmi,MF,MC,L,G,B))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
			if (MAXHITS <=Hits || MAX_MISMATCHES == 0) continue;
		}


		if(Hits=Map_Strand_Guess(MAX_MISMATCHES,MAXHITS,L,fwfmi,revfmi,G.Guessed))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guessed,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
		}
		else if(Hits=Map_Strand_Guess(MAX_MISMATCHES,MAXHITS,L,fwfmi,revfmi,G.Guess_Complement))
		{
			Mapped++;//=Hits;
			//Print_Hits(G.Guess_Complement,L,Output_File,ONEFMINDEX,revfmi,fwfmi,O);
			Print_Hits_To_BAT(G.Guess_Complement,L,Output_File,ONEFMINDEX,revfmi,F.FILETYPE);
		}
		else
		{
			//printf ("%u\t%s\n",Actual_Tag,R.Tag_Copy);
		}

	}
	char End_Mark='&';
	fwriteX(&End_Mark,1,1,Output_File);//write sentinel..
	return Mapped;
}*/

//{-----------------------------------  MAPPING MODES ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Map_Strand
 *  Description:  Do the mapping, along the given strand.. 
 * =====================================================================================
 */
int Map_Strand(int & Last_Mis,bool batmeth1,int MAX_MISMATCHES,int MAXHITS,LEN & L,BWT *fwfmi,BWT* revfmi,MEMX & M)
{
	unsigned Hits=0,H;
	char* Current_Tag=M.Current_Tag;
	M.Strand=Current_Tag[L.STRINGLENGTH];
	M.Larger_Than_Ten=0;
	if(H=Zero_Mismatch(batmeth1,Current_Tag,L,revfmi,M,fwfmi))
	{
		Hits+=H;Last_Mis=0;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==0 ) return Hits;
	if (H=One_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;Last_Mis=1;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==1 ) return Hits;
	if (H= Two_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;Last_Mis=2;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==2 ) return Hits;
	if (H=Three_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;Last_Mis=3;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==3 ) return Hits;
	if (H= Four_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;Last_Mis=4;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==4 ) return Hits;
	if (H= Five_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;Last_Mis=5;if (MAXHITS <= Hits) return Hits;
	}
	return Hits;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Map_Strand_Guess
 *  Description:  Do the mapping, assuming the strand was guessed..
 * =====================================================================================
 */
unsigned Map_Strand_Guess(bool batmeth1,int MAX_MISMATCHES,int MAXHITS,LEN & L,BWT *fwfmi,BWT* revfmi,MEMX & M)
{
	unsigned Hits=0,H;
	char* Current_Tag=M.Current_Tag;
	M.Strand=Current_Tag[L.STRINGLENGTH];
	M.Larger_Than_Ten=0;
	if(!M.Guessed && (H=Zero_Mismatch(batmeth1,Current_Tag,L,revfmi,M,fwfmi)))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	else
	{
		M.Left_Mishits_Pointer=0;M.Right_Mishits_Pointer=0;M.Possible_20_Pointer=0;M.Possible_02_Pointer=0;M.Mismatches_Forward_Pointer=0;M.Mismatches_Backward_Pointer=0;M.Two_Mismatches_At_End_Pointer=0;M.Two_Mismatches_At_End_Forward_Pointer=0;
	}
	if (MAX_MISMATCHES ==0 ) return Hits;
	if (H=One_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==1) return Hits;
	if (H= Two_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==2 ) return Hits;
	if (H=Three_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==3 ) return Hits;
	if (H= Four_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	if (MAX_MISMATCHES ==4 ) return Hits;
	if (H= Five_Mismatch(batmeth1,Current_Tag,L,MAXHITS,fwfmi,revfmi,M))
	{
		Hits+=H;if (MAXHITS <= Hits) return Hits;
	}
	return Hits;
}
//}-----------------------------------  MAPPING MODES ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Split_File
 *  Description:  split file into C pieces..
 * =====================================================================================
 */
/*off64_t Split_File(FILE* F, int C, char FILETYPE)
{
	char Description[MAXDES];
	off64_t File_Size,Last,Split;

	off64_t Current=ftello64(F);//mark last comment...
	fseek(F, 0L, SEEK_END);
off64_t X=ftello64(F);
	Split = ftello64(F)/C;
	fseek(F, Split, SEEK_SET);
	if (FILETYPE == FQ)
	{
		for(;;)//ignore comments...
		{
			//gzgets(Input_File,Description,MAXDES);
			fgets(Description,MAXDES,F);
			Last=ftello64(F);//mark last comment...
			if (Description[0] != '@') break;
		}
	}
	fseek(F,Current,SEEK_SET);//go to first reads..
	return Last;
}*/
