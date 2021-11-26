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

const int BUFSIZE=500000000;
FILE* File_Open(const char* File_Name,const char* Mode);
char Fill_Char;
int main(int argc, char* argv[])
{
	int Buf_Size=BUFSIZE;
	string BaseName;
	char Rand_Array[]={'a','c','g','t'};
	
	srand(0);
	if (argc>1)
	{
		FILE* INFILE=File_Open(argv[1],"r");
		BaseName=argv[1];
		string N_File=BaseName+".len";
		string Bin_File=BaseName+".bin";
		FILE* NFILE=File_Open(N_File.c_str(),"w");
		FILE* BIN=File_Open(Bin_File.c_str(),"w");
		unsigned Offset=0;

		char Convert[255];
        for (int i=0;i<255;i++) {Convert[i]='c';}
        Convert['a']=Convert['A']='a';
        Convert['g']=Convert['G']='g';
        Convert['t']=Convert['T']='t';
        Convert['\r']='\r';Convert['\n']='\n';Convert['-']='-';Convert['N']=Convert['n']='-';

		if (argc>2) Buf_Size=atoi(argv[2]);
		char* Buffer=new char[Buf_Size];
		Buffer[Buf_Size]=0;

		int i;
		while (fgets(Buffer,Buf_Size,INFILE))
		{
			
			if (Buffer[0] != '>') 
			{
				for(i=0;Buffer[i] && Buffer[i]!='\n' && Buffer[i]!='\r' && Buffer[i]!='\t' && Buffer[i]!=' ' && i<Buf_Size ;i++)
				{
					Buffer[i]=Convert[Buffer[i]];
					fputc((int)Buffer[i],BIN);
					Offset++;	
				}
			}
			else
			{
                                for(i=0;Buffer[i] && Buffer[i]!='\n' && Buffer[i]!='\r' && Buffer[i]!='\t' && Buffer[i]!=' ' && i<Buf_Size ;i++);
                                Buffer[i] = '\0';
                                if(Offset == 0) {
                                    fprintf(NFILE,"%s",Buffer+1);
                                }
				else fprintf(NFILE,"\t%u\n%s",Offset, Buffer+1);
				Offset=0;
			}
			//printf("%s",Buffer);
		}
                fprintf(NFILE,"\t%u\n",Offset);
	}
	else
	{
		printf ("Filter file out of ambiguos bases..\n");
		printf ("Command Line: Filter file_name..\n");
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
