#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <string>
#include <string.h>

using namespace std;

struct N_Int
{
	unsigned Start;
	unsigned End;
};

const int BUFSIZE=60000;
FILE* File_Open(const char* File_Name,const char* Mode);
int main(int argc, char* argv[])
{
	int Buf_Size=BUFSIZE;
	bool In_Amb=false;
	bool RANDMODE=false;
    char seps[] = " \t\n\r";
    char* token;
	string genomefa = "";
    string prefix = "";
    string HELP = "genome2cg -g genome.fa\n" \
                  " -g genome.fa\n" \
                  " -p prefix\n" \
                  " -b buffer size for max line\n";
	char Rand_Array[]={'a','c','g','t'};
    for(int i=1; i<argc; i++){
        if(strcmp(argv[i], "-g")==0){
            genomefa = argv[++i];
        }else if(strcmp(argv[i], "-p")==0){
            prefix = argv[++i];
        }else if(strcmp(argv[i], "-b")==0){
            Buf_Size = atoi(argv[++i]);
        }else{
            fprintf(stderr, "%s\n", HELP.c_str());
        }
    }
    if(prefix == ""){
        prefix = genomefa;
    }
	
	srand(0);
	if (argc>2)
	{
		FILE* INFILE=File_Open(genomefa.c_str(),"r");
        string Genome_File = prefix+".batmeth2.fa";
        FILE* GFILE = File_Open(Genome_File.c_str(), "w");

		unsigned Offset=0;
		char ConvertCT[255]; char ConvertGA[255];

        {
            fprintf(stderr,"Converting C->T\n");
            for (int i=0;i<255;i++) {ConvertCT[i]='t';}
            ConvertCT['C']=ConvertCT['c']='T';ConvertCT['a']=ConvertCT['A']='A';ConvertCT['g']=ConvertCT['G']='G';
            ConvertCT['t']=ConvertCT['T']='T';ConvertCT['\r']='\r';ConvertCT['\n']='\n';
            ConvertCT['-']='-';
            ConvertCT['N']=ConvertCT['n']='-';
        }//C-->T
        {
            fprintf(stderr,"Converting G->A\n");
            for (int i=0;i<255;i++) {ConvertGA[i]='a';}
            ConvertGA['a']=ConvertGA['A']='A';ConvertGA['G']=ConvertGA['g']='A';
            ConvertGA['t']=ConvertGA['T']='T';ConvertGA['C']=ConvertGA['c']='C';
            ConvertGA['\r']='\r';ConvertGA['\n']='\n';ConvertGA['-']='-';
            ConvertGA['N']=ConvertGA['n']='-';
        }// GA

		char* Buffer=new char[Buf_Size];
		Buffer[Buf_Size]=0;
        char* chrom = new char[Buf_Size];

        long genomeFoffset = ftell(INFILE);
        int process = -1; //0 CT, 1 GA

		int i;
		while (fgets(Buffer,Buf_Size,INFILE))
		{
			
			if (Buffer[0] != '>') 
			{
				for(i=0;Buffer[i] && i<Buf_Size ;i++) // Buffer[i]!='\n' && Buffer[i]!='\r'
				{
                    if(process == 0) Buffer[i]=ConvertCT[Buffer[i]];
                    else Buffer[i]=ConvertGA[Buffer[i]];
					fputc((int)Buffer[i],GFILE);
				}
			}
			else
			{
                if(process == -1){ // C2T
                    token = strtok(Buffer+1, seps);
				    fprintf(GFILE,">f%s\n",token);
                    genomeFoffset = ftell(INFILE); // remeber offset, 此时记载的是行尾
                    strcpy(chrom, token);
                    process = 0;
                    continue;
                }
                if(process == 0){ // G2A
                    fprintf(GFILE,">r%s\n",chrom);
                    fseek(INFILE, genomeFoffset, SEEK_SET);
                    process = -1;
                }
			}
		}
        //last chrom, g2a
        fprintf(GFILE,">r%s\n",chrom);
        fseek(INFILE, genomeFoffset, SEEK_SET);
        process = -1;
        while(fgets(Buffer,Buf_Size,INFILE)){
            for(i=0;Buffer[i] && i<Buf_Size ;i++) // Buffer[i]!='\n' && Buffer[i]!='\r'
            {
                if(process == 0) Buffer[i]=ConvertCT[Buffer[i]];
                else Buffer[i]=ConvertGA[Buffer[i]];
                fputc((int)Buffer[i],GFILE);
            }
        }
        fclose(INFILE);
        fclose(GFILE);
	}
	else
	{
		fprintf(stderr, "%s\n", HELP.c_str());
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
