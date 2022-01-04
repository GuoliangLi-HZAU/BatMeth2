#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <map>
#define CHROMSIZE 100
#define BATBUF 2000
/*
	DMC count annotation in TE/LTR/LINE/SINE/GENE/EXON/5'UTR/3'UTR ...
*/
struct GenomeDMC
{
	int *plusDMC;
	int *NegDMC;
	int Index;
};
struct GenomeAllC
{
//	char* Genome;
	int *plusAllC;
	int *NegAllC;
	int Index;
};
using namespace std;
map <string,int> String_Hash;
void Show_Progress(float Percentage);
void Splitpath(string & str, string & fname);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....

FILE* File_Open(const char* File_Name,const char* Mode);
unsigned Total_Reads;
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
unsigned totalC=0;
unsigned totalDMC=0;
int main(int argc,char* argv[]){
	//Help
	time_t Start_Time,End_Time;
	printf("\nBatMeth2:  DMCannotationPlot v1.0\n");
	const char* Help_String="Command Format :   DMCannotationPlot [options] -o <Out_File> -G GENOME -g <GFF files..> -d <dmc file> -c<mC context>\n"
		"\nUsage:\n"
		"\t-o|--out       Output file name.\n"
		"\t-G|--genome    Genome file\n"
		"\t-d|--dmcFile   dmc file. Format: Chrome Location strand\n"
		"                 Format: chr pos strand \n"
		"\t-g|--gff       Gff files, 1 or more\n"
		"\t-c|--context   mC context, CG[default]/CHG/CHH/C. \n"
		"\t-h|--help";
	
	int GffInFileStart=0;
	int GffInFileEnd=0;
	char *dmcFile;
	string Geno;
	int Genome_CountX=0;
	char* Output_Name;
	string context="CG";
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-d") ||!strcmp(argv[i], "--dmcFile")  )
		{
			dmcFile=argv[++i];
		}else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--gffFile"))
		{
			GffInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {GffInFileEnd=--i;}else {GffInFileEnd=i ;}
		}else if(!strcmp(argv[i], "-G") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile"))
		{
			Output_Name=argv[++i];
		}else if(!strcmp(argv[i], "-c") || !strcmp(argv[i], "--context"))
		{
			context=argv[++i];
		}else
		{
			printf("%s \n",Help_String);
			exit(0);
		}
		
	}
	if(argc<=1) 
	{	
		printf("%s \n",Help_String);
		exit(0);
	}
	try
	{
		time(&Start_Time);
		// Load Genome and initial Hash
		string G=Geno;G+=".bin";
		FILE* BINFILE=File_Open(G.c_str(),"r");
		
		string L=Geno;L+=".ann.location";
		FILE* Location_File=File_Open(L.c_str(),"r");
		
		fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
		char* Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n"); 
		if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");
			
		printf("LocationFile: %s\n",L.c_str() );
		struct Offset_Record
		{
			char Genome[40];
			unsigned Offset;
		} Temp_OR; 
	
		while (fgets(Temp_OR.Genome,39,Location_File)!=0)//count genomes..
		{
			fgets(Temp_OR.Genome,39,Location_File);
			Genome_CountX++;	
		}
		rewind(Location_File);
		Offset_Record Genome_Offsets[Genome_CountX+1];
		
		int Genome_Count=0;
		while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<=Genome_CountX)
		{
			Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
			fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
			for(int i=0;i<40;i++) 
			{
				if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
				{ 
					Genome_Offsets[Genome_Count].Genome[i]=0;
					break;
				} 
			}
			Genome_Count++;	
		}
		Genome_Count--;
		
		char* Split_Point=Org_Genome;//split and write...
		GenomeDMC Methy_DMC[Genome_Count];
		GenomeAllC Methy_AllC[Genome_Count];
		for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
		{
			printf("%s, ",Genome_Offsets[i].Genome);
			String_Hash[Genome_Offsets[i].Genome]=i;
			//meth init
			Methy_DMC[i].plusDMC =new int[Genome_Offsets[i+1].Offset];
			Methy_DMC[i].NegDMC =new int[Genome_Offsets[i+1].Offset];
			Methy_DMC[i].Index=i;
			//allC init
			//Methy_AllC[i].Genome=(Split_Point+=Genome_Offsets[i].Offset);
			Methy_AllC[i].plusAllC = new int[Genome_Offsets[i+1].Offset];
			//Methy_AllC[i].NegAllC = new int[Genome_Offsets[i+1].Offset];
			Methy_AllC[i].Index=i;
			//
			char* Genome=(Split_Point+=Genome_Offsets[i].Offset); //Methy_AllC[i].Genome;
			for(int l=0;l<Genome_Offsets[i+1].Offset;l++)
			{
				char genome_Char = toupper(Genome[l]);
				char genome_CharFor1,genome_CharFor2;
				if(l+2<Genome_Offsets[i+1].Offset)
				{
					 genome_CharFor1= toupper(Genome[l+1]);
					 genome_CharFor2 = toupper(Genome[l+2]);
				}
				char genome_CharBac1,genome_CharBac2;
				if(l> 2)
				{
					genome_CharBac1 = toupper(Genome[l-1]);
					genome_CharBac2 = toupper(Genome[l-2]);
				}
				if(genome_Char=='C' /*|| genome_Char=='G'*/ ) totalC++;
				if(context=="CG")
				{
					if(genome_Char=='C' && genome_CharFor1== 'G') Methy_AllC[i].plusAllC[l]=1;
					//if(genome_Char=='G' && genome_CharBac1== 'C') Methy_AllC[i].NegAllC[l]=1;
				}else if(context=="CHG")
				{
					if(genome_Char=='C' && genome_CharFor1!= 'G' && genome_CharFor2=='G' ) Methy_AllC[i].plusAllC[l]=1;
					//if(genome_Char=='G' && genome_CharBac1!= 'C' && genome_CharBac2=='C' ) Methy_AllC[i].NegAllC[l]=1;
				}else if(context=="CHH")
				{
					if(genome_Char=='C' && genome_CharFor1!= 'G' && genome_CharFor2!='G' ) Methy_AllC[i].plusAllC[l]=1;
					//if(genome_Char=='G' && genome_CharBac1!= 'C' && genome_CharBac2!='C' ) Methy_AllC[i].NegAllC[l]=1;
				}else if(context=="C")
				{
					if(genome_Char=='C') Methy_AllC[i].plusAllC[l]=1;
					//if(genome_Char=='G') Methy_AllC[i].NegAllC[l]=1;
				}else
				{
					fprintf(stderr,"\nPlease assign the mC context used argument -c|--context.\n");
					exit(0);
				}
				
			}
			
		}
		printf("Loaded\n");
		
		//load dmc file
		printf("\nLoading dmc file...\n");
		FILE* DMCFILE=File_Open(dmcFile,"r");
		int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
		fseek(DMCFILE, 0L, SEEK_END);off64_t File_Size=ftello64(DMCFILE);rewind(DMCFILE);
		Init_Progress();
		char DMCbuffer[BATBUF],Chrom[CHROMSIZE],Strand;unsigned pos;
		
		while (!feof(DMCFILE)) 
		{
			Total_Reads++;
			Progress++;
			if (Progress>Number_of_Tags) 
			{
				off64_t Current_Pos=ftello64(DMCFILE);
				off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
				Number_of_Tags=(File_Size/Average_Length)/20;
				Progress=0;
				Show_Progress(Current_Pos*100/File_Size);
			}
			fgets(DMCbuffer,BATBUF,DMCFILE);
			sscanf(DMCbuffer,"%s%u%s",Chrom,&pos,&Strand);
			if(String_Hash.count(Chrom)>0)
			{
				int H=String_Hash[Chrom];
				pos--;//0-- for array start 0
				if(!(pos>=0 && pos < Genome_Offsets[H+1].Offset )) 
				{
					printf("\nout of limit of chrome length, Line : %d\n",Progress);
					exit(0);
				}
				totalDMC++;
				if(Strand=='+')
				{
					Methy_DMC[H].plusDMC[pos]=1;
				}else if(Strand=='-')
				{
					Methy_DMC[H].NegDMC[pos]=1;
				}else printf("wrong strand...\n");
			}
		}
		Done_Progress();
		fclose(DMCFILE);
		printf("Methratio file Loaded.\n");
		
		FILE* OUTFILE=File_Open(Output_Name,"w"); // store OutPut
		//load Gff file
		for(int f=GffInFileStart;f<=GffInFileEnd;f++){
			printf("\nProcessing %d out of %d. GffFile: %s\n", f-GffInFileStart+1,GffInFileEnd-GffInFileStart+1, argv[f]);
			FILE* GFFINFILE=File_Open(argv[f],"r");
			Progress=0;Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(GFFINFILE, 0L, SEEK_END);off64_t File_Size=ftello64(GFFINFILE);rewind(GFFINFILE);//
			Init_Progress();Total_Reads=0;
			char Gff[BATBUF],temp[BATBUF];
			string Symbol;
			unsigned start=0,end=0;
			unsigned DMCcount=0;
			unsigned AllCcount=0;
			while(!feof(GFFINFILE))
			{
				Total_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(GFFINFILE);
					off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
					Number_of_Tags=(File_Size/Average_Length)/20;
					Progress=0;
					Show_Progress(Current_Pos*100/File_Size);
				}
				fgets(Gff,BATBUF,GFFINFILE);
				sscanf(Gff,"%s%s%s%u%u%s%s%s%s",Chrom,temp,temp,&start,&end,temp,&Strand,temp,Symbol.c_str());
				
				if(String_Hash.count(Chrom)>0)
				{
					int H=String_Hash[Chrom];
					if(start>0 && end <= Genome_Offsets[H+1].Offset )
					{
						for(int i=start-1;i<end;i++) //-1 for array start 0
						{
							//if(Strand=='+' || Strand=='.') 
							{
								DMCcount+=Methy_DMC[H].plusDMC[i];
								AllCcount+=Methy_AllC[H].plusAllC[i];
							}
							/*
							if(Strand=='-' || Strand=='.' ) 
							{
								DMCcount+=Methy_DMC[H].NegDMC[i];
								AllCcount+=Methy_AllC[H].NegAllC[i];
							}
							*/
						}
					}
				}
				
			}
			Done_Progress();
			fclose(GFFINFILE);
			string str= argv[f];string fname;
			Splitpath(str, fname);
			fprintf(OUTFILE,"%s\t%d\tDMC\n",fname.c_str(),DMCcount);
			fprintf(OUTFILE,"%s\t%d\tAllC\n",fname.c_str(),AllCcount);
		}
		fprintf(OUTFILE,"Genome\t%d\ttotalC\n",totalC);
		fprintf(OUTFILE,"Genome\t%d\ttotalDMC\n",totalDMC);
		time(&End_Time);
		fclose(OUTFILE);
		printf("\nTime Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
		exit(0);
	}	
	catch(char* Err)
	{
		printf(Err);
		exit(-1);
	}
	return 0;
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

void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}

void Splitpath(string & str, string & fname)
{
	int m = str.rfind("."); //reverse find  ��.�� ,return position(>=0) of '.' 
	int n = str.rfind("\\");
	if(n==string::npos ) n = str.rfind("/");
	//string dir = str.substr(0,n); //lacation of the file //substr��(start,length) 
	//string name = str.substr(n+1,str.length()-n-1); //file name
	if(n==string::npos ) fname= str.substr(0,m-n-1);
	else if(m!=string::npos) fname= str.substr(n+1,m-n-1); //file name remove suffix 
	else fname=str;
}
