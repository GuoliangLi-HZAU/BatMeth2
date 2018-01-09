/*
 * ReverseComplteFQ.cpp
 *
 *  Created on: 2015/12/02
 *  Author: biozc
 */
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXTAG 500
#define FQ       2
#define FA       3
#define TRUE    1
#define FALSE    0
struct INFILE
{
	FILE* Input_File;
	char FILETYPE;
};
struct READ
{
	char Description[MAXTAG];
	char Tag_Copy[MAXTAG];
	char Quality[MAXTAG];
	unsigned Read_Number;
	int Real_Len;
	char Plus[MAXTAG];
};
char Char_To_CharC[256];
FILE* File_Open(const char* File_Name,const char* Mode);
void Analyze_File(INFILE & I);
char Read_Tag(FILE *Input_File,const char FILETYPE, READ & Read );
void Reverse_Quality(char* Dest,const READ & R,int StringLength);

int main(int argc, char* argv[]) {

	if(argc<3)
	{
		printf("Usage: ReverseComplteFQ Input.fq Out.fq\n");
		exit(0);
	}
	INFILE Head_File;
	Head_File.Input_File=File_Open(argv[1],"r");
	READ R;READ outR;
	Analyze_File(Head_File);
	//--reverce
	Char_To_CharC['A']='T';Char_To_CharC['C']='G';
	Char_To_CharC['G']='C';Char_To_CharC['T']='A';
	Char_To_CharC['a']='t';Char_To_CharC['c']='g';
	Char_To_CharC['g']='c';Char_To_CharC['t']='a';
	Char_To_CharC['n']='n';Char_To_CharC['N']='N';
	//OutFile
	FILE* OutFile=File_Open(argv[2],"w");
	while (Read_Tag(Head_File.Input_File,Head_File.FILETYPE,R))
	{
		R.Real_Len=0;
		for(;R.Tag_Copy[R.Real_Len]!=0 && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
		if(R.Tag_Copy[R.Real_Len] == '\n') R.Tag_Copy[R.Real_Len]='\0';
		if(R.Quality[R.Real_Len] == '\n') R.Quality[R.Real_Len]='\0';
		outR=R;
		int i;
		for (i=0;i<=R.Real_Len-1;i++){outR.Tag_Copy[R.Real_Len-1-i]=Char_To_CharC[R.Tag_Copy[i]];};
		Reverse_Quality(outR.Quality,R,R.Real_Len);
		fprintf(OutFile,"%s\n%s\n+\n%s\n",outR.Description,outR.Tag_Copy,outR.Quality);
	}
	fclose(Head_File.Input_File);
	fclose(OutFile);
}

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
	Handle=fopen(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....\n",File_Name);
		exit(-1);
	}
	else return Handle;
}

void Analyze_File(INFILE & I)
{
	char Description[MAXTAG];
	char Tag_Copy[MAXTAG];
	char Current_Tag[MAXTAG];
	char Quality[MAXTAG];
	char Plus[MAXTAG];
	unsigned Last=0;

	fseek(I.Input_File, 0L, SEEK_END);
	rewind(I.Input_File);

	for(;;)//ignore comments...
	{
		if(!fgets(Description,MAXTAG,I.Input_File)) {printf("Analyze_File(): error reading file...\n");exit(-1);}
		if (Description[0] != '#') break;
		Last=ftello(I.Input_File);//mark last comment...
	}
	if(!fgets(Current_Tag,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}

	if(!fgets(Quality,MAXTAG,I.Input_File)){printf("Analyze_File(): error reading file...\n");exit(-1);}//plus
	if (Quality[0]=='>') I.FILETYPE=FA;
	else
	{
		I.FILETYPE=FQ;
		if(Quality[0] != '+' || Description[0] != '@') {printf("Init_Variables: Cannot determine file type ...\n");exit(1);}
	}
	fseek(I.Input_File,Last,SEEK_SET);//go top
}
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
	if (fgets(Read.Description,MAXTAG,Input_File)!=0)// read a tag...
	{
		char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
		//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
		if(!fgets(Read.Tag_Copy,MAXTAG,Input_File)) {printf ("Read_Tag():Error reading file Tag..\n");exit(-1);};//tag
		if (FILETYPE == FQ)
		{
			//gzgets(Input_File,Plus,MAXTAG);//plus
			if(!fgets(Read.Plus,MAXTAG,Input_File)){ printf ("Read_Tag():Error reading file Plus..\n");exit(-1);};//plus
			//gzgets(Input_File,Quality,MAXTAG);//phred
			if(!fgets(Read.Quality,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file Qulity..\n");exit(-1);};//phred
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

void Reverse_Quality(char* Dest,const READ & R,int StringLength)
{
	const char* Quality =R.Quality;
	for (int i=StringLength-1;i>=0;i--)
	{
		*Dest=Quality[i];Dest++;
	}
	*Dest=0;
}
