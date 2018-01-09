#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <string>

using namespace std;

struct N_Int
{
	unsigned Start;
	unsigned End;
};
int chromname_len=400;
const int BUFSIZE=60000;
FILE* File_Open(const char* File_Name,const char* Mode);
char Fill_Char;
int main(int argc, char* argv[])
{
	int Buf_Size=BUFSIZE;
	bool In_Amb=false;
	bool RANDMODE=false;
	string BaseName;
	char Rand_Array[]={'a','c','g','t'};
	
	srand(0);
	if (argc>1)
	{
		FILE* INFILE=File_Open(argv[1],"r");
		BaseName=argv[1];

		char* Buffer=new char[Buf_Size];
		Buffer[Buf_Size]=0;

		while (fgets(Buffer,Buf_Size,INFILE))
		{
			
			if (Buffer[0] != '>') 
			{
				printf("%s",Buffer);
			}
			else if (Buffer[0] == '>') 
			{
				char tmpstring[chromname_len];
				sscanf(Buffer,"%[^ |^\t|^\n]",tmpstring);
				printf("%s\n",tmpstring);
			}
			
		}
	}
	else
	{
		printf ("Command Line: genome_filter file_name > outgenome..\n");
	}
	//for (int i=0;i<NList_Index;i++) fwrite(&NList[i],1,sizeof(N_Int),NFILE);

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
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}
