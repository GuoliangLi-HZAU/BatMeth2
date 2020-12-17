#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include "limits.h"
#include <map>
#include <sstream>
#define CHROMSIZE 1000
#define BATBUF 10000

struct Methy_Hash
{
	int *plusCount_C,*plusCount_CT;
	int *NegCount_C,*NegCount_CT;
	char *plusMethContext;
	char *NegMethContext;
	int Index;
//	int *binsPlusCount_C,*binsPlusCount_CT;
//	int *binsNegCount_C,*binsNegCount_CT;
};
struct GeneDensity
{
	unsigned *PN_Cover;
	unsigned *plusCover;
	unsigned *NegCover;
	int Index;
};
struct Methy_Gff
{
	long *CG_C,*CHG_C,*CHH_C;
	long *CG_CT,*CHG_CT,*CHH_CT; // density plot
	long AverPerCG,AverPerCHG,AverPerCHH; //DM
	long AverCG,AverCHG,AverCHH;
	int countCG,countCHG,countCHH;
		//int *MethContext;
	int Index;
};
using namespace std;
bool Collision=false;
map <string,int> String_Hash;
map <string,int> Context_Hash;


//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List);
void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List,GeneDensity & GeneD_List);
void caculateHeatmap(const char* type,int start,int end,Methy_Hash MethyList,char Strand,const char* id,FILE *methGffcg,FILE *methGffchg,FILE *methGffchh,char* chrom, int geneS,int geneE,bool printtitle=true); //FILE *methGffcg_matrix,FILE *methGffchg_matrix,FILE *methGffchh_matrix ,
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....
unsigned u=0;
unsigned Total_Reads;//total hits in the input file...
float binspan=0.01;
float bodyspan = 0.01;
unsigned nLevel=ceil(1/(double)binspan)-1;//
unsigned bodynLevel = nLevel;
unsigned bodyLen=2000;
unsigned longerlen = nLevel;
string int2str(int &int_temp)
{
        stringstream stream;
        stream<<int_temp;
        return stream.str();
}
string unsigned2str(unsigned &int_temp)
{
        stringstream stream;
        stream<<int_temp;
        return stream.str();   
}
void str2int(int &int_temp,string &string_temp)
{
	stringstream stream(string_temp);
	stream>>int_temp;
}
//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
long AverPerCG=0,AverPerCHG=0,AverPerCHH=0;
long AverCG=0,AverCHG=0,AverCHH=0;
int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	printf("\nBatMeth2:  MethyGff v1.0\n");
	const char* Help_String="Command Format :   methyGff [options] -o <OUT_PREFIX> -G GENOME -gff <GFF file>/-gtf <GTF file>/-b <bed file>/-b4 <bed4 file> -m <from Split methratio outfile> [-B] [-P]\n"
		"\nUsage:\n"
		"\t-o|--out         Output file prefix\n"
		"\t-G|--genome      Genome\n"
		"\t-m|--methratio   Methratio output file.\n"
		"\t-c|--coverage    >= <INT> coverage. default:4\n"
		"\t-C               <= <INT> coverage. default 600.\n"
		"\t-nC              >= <INT> Cs per bins or genes. default:1\n"
		"\t-gtf|-gff        Gtf/gff file\n"
		"\t-b|--BED         Bed file, chrom start end\n"
		"\t-b4              Bed file, chrom start end strand\n"
		"\t-b5              Bed file, chrom start end geneid strand\n"
		//"\t--bins           DMR bins step, default 1000bp.\n"
		"\t-d|--distance    DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000\n"
		"\t-B|--body        For different analysis input format, gene/TEs body methylation level. [Different Methylation Gene(DMG/DMT...)]\n"
		"\t-P|--promoter    For different analysis input format.[Different Methylation Promoter(DMP)]\n"
		"\t--TSS            Caculate heatmap for TSS. [Outfile: outPrefix.TSS.cg.n.txt]\n"
		"\t--TTS            Caculate heatmap for TTS. [Outfile: outPrefix.TTS.cg.n.txt] \n"
		"\t--GENE           Caculate heatmap for GENE and flank 2k. [Outfile: outPrefix.GENE.cg.n.txt] \n"
		"\t-s|--step        Gene body and their flanking sequences using an overlapping sliding window of 2% of the sequence length at a step of 1% of the sequence length. So default step: 0.01 (1%)\n"
		"\t-bl|--bodyLen    Body length to which all regions will be fit. (default: same as -d)\n"
		"\t-S|--chromStep   Caculate the density of genes/TEs in chromsome using an overlapping sliding window of 100000bp at a step of 50000bp, must equal \"-s\" in Split.. default step: 50000(bp)\n"//
		"\t-h|--help";

	int Genome_CountX=0;
	char* Output_Name;
	//char* methInfileName;
	int InFileStart=0,InFileEnd=0;
	string Geno;
	char *InFile;
	bool InputGff=false;
	int GffInFileStart=0,GffInFileEnd=0;
	bool InputBed=false;
	int BedInFileStart=0,BedInFileEnd=0;
	int PrefixStart=0,PrefixEnd=0;
	int distance=2000;
	int distanceHeatmap=2000;
	int chromStep=50000;
	int binsStep=1000;
	int maxcover=1000;
	bool Diff=false;
	bool PU=false;
	bool TSS=false;
	bool TTS=false;
	bool GENE=false;
	bool GTF=false;
	int Coverage = 4;
	int nC=1;
	int beddistance=1;
	bool bed4 = false, bed5=false; int showbed = 0;
	//
	if(argc<8)
	{
		printf("%s \n",Help_String);
		exit(0);
	}
	for(int i=1;i<argc;i++)
	{
		if(!strcmp(argv[i], "-o") ||!strcmp(argv[i], "--out")  )
		{
			Output_Name=argv[++i];
			PrefixStart=i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {PrefixEnd=--i;}else {PrefixEnd=i ;}
		}
		else if(!strcmp(argv[i], "-G") || !strcmp(argv[i], "--genome"))
		{
			Geno=argv[++i];
		}else if(!strcmp(argv[i], "-c") || !strcmp(argv[i], "--coverage"))
		{
			Coverage=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-C")){
			maxcover=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-nC"))
		{
			nC=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-s") || !strcmp(argv[i], "--step"))
		{
			binspan=atof(argv[++i]);
			nLevel=ceil(1/(double)binspan)-1;
			printf("binspn: %f, nLevel: %d\n", binspan, nLevel);
		}else if(!strcmp(argv[i], "-bl") || !strcmp(argv[i], "--bodyLen"))
		{
			bodyLen=atoi(argv[++i]);
			bodyspan = binspan/((double)bodyLen/distance);
			bodynLevel=ceil(1/(double)bodyspan)-1;
		}
		else if(!strcmp(argv[i], "--bins"))
		{
			binsStep=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-B") || !strcmp(argv[i], "--body"))
		{
			Diff=true;
		}
		else if(!strcmp(argv[i], "-P") || !strcmp(argv[i], "--promoter"))
		{
			PU=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--TSS"))
		{
			TSS=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--TTS"))
		{
			TTS=true;//promoter UP
		}
		else if( !strcmp(argv[i], "--GENE"))
                {
                        GENE=true;//promoter UP
                }
		else if(!strcmp(argv[i], "-S") || !strcmp(argv[i], "--chromStep"))
		{
			chromStep=atoi(argv[++i]);
		}
		else if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance"))
		{
			distance=atoi(argv[++i]);
			bodyLen = distance;
			distanceHeatmap = distance;
		}
		else if(!strcmp(argv[i],"-m") || !strcmp(argv[i],"--methratio"))
		{
			InFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {InFileEnd=--i;}else {InFileEnd=i ;}
			//methInfileName=argv[++i];
		}
		else if(!strcmp(argv[i], "-gtf") || !strcmp(argv[i], "-gff"))
		{
			InputGff=true;
			if(!strcmp(argv[i], "-gtf")) GTF=true;
			GffInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {GffInFileEnd=--i;}else {GffInFileEnd=i ;}
		}else if(!strcmp(argv[i], "-b") || !strcmp(argv[i], "--BED"))
		{
			InputBed=true;
			BedInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {BedInFileEnd=--i;}else {BedInFileEnd=i ;}
		}else if(!strcmp(argv[i], "-b4") || !strcmp(argv[i], "--BED4"))
		{
			InputBed=true;
			BedInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {BedInFileEnd=--i;}else {BedInFileEnd=i ;}
			showbed = 1;
			bed4 = true;
		}else if(!strcmp(argv[i], "-b5") || !strcmp(argv[i], "--BED5"))
		{
			InputBed=true;
			BedInFileStart=++i;
			while(i!=(argc-1) && argv[i][0]!='-'){i++;continue;}
			if(argv[i][0]=='-') {BedInFileEnd=--i;}else {BedInFileEnd=i ;}
			bed5 = true; showbed = 1;
		}
		else
		{
			printf("%s \n",Help_String);
			exit(0);
		}
	}
	if (argc >1) 
	{
		try
		{
			if(bodyLen != distance){
				bodyspan = binspan/((double)bodyLen/distance);
				bodynLevel=ceil(1/(double)bodyspan)-1;
			}else{
				bodyspan=binspan;
				bodynLevel = nLevel;
			}

			printf("Output prefix: %s\n", Output_Name);
			time(&Start_Time);

			string L=Geno;L+=".ann.location";

			FILE* Location_File=File_Open(L.c_str(),"r");
			
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

			GeneDensity GeneD_List[Genome_Count];
			Methy_Hash Methy_List[Genome_Count];
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				String_Hash[Genome_Offsets[i].Genome]=i;
				//meth ini
				Methy_List[i].plusMethContext =new char[Genome_Offsets[i+1].Offset]();
				Methy_List[i].NegMethContext =new char[Genome_Offsets[i+1].Offset]();
				Methy_List[i].plusCount_C = new int[Genome_Offsets[i+1].Offset]();
				Methy_List[i].plusCount_CT = new int[Genome_Offsets[i+1].Offset]();
				Methy_List[i].NegCount_C = new int[Genome_Offsets[i+1].Offset]();
				Methy_List[i].NegCount_CT = new int[Genome_Offsets[i+1].Offset]();
				Methy_List[i].Index=i;
				//Methy_List[i].binsPlusCount_C =new int[(int)ceil((double)Genome_Offsets[i+1].Offset/binsStep)];
				//Methy_List[i].binsPlusCount_CT =new int[(int)ceil((double)Genome_Offsets[i+1].Offset/binsStep)];
				//Methy_List[i].binsNegCount_C =new int[(int)ceil((double)Genome_Offsets[i+1].Offset/binsStep)];
				//Methy_List[i].binsNegCount_CT =new int[(int)ceil((double)Genome_Offsets[i+1].Offset/binsStep)];
				
				unsigned len=ceil(Genome_Offsets[i+1].Offset/double(chromStep))-1;
				GeneD_List[i].PN_Cover=new unsigned[Genome_Offsets[i+1].Offset]();
				GeneD_List[i].plusCover=new unsigned[Genome_Offsets[i+1].Offset]();
				GeneD_List[i].NegCover=new unsigned[Genome_Offsets[i+1].Offset]();
				GeneD_List[i].Index=i;
			}
			printf("Loaded\n");
			string CG="CG",CHG="CHG",CHH="CHH";
			Context_Hash[CG.c_str()]=1;Context_Hash[CHG.c_str()]=2;Context_Hash[CHH.c_str()]=3;
			
			unsigned pos;float methratio=0;int countC=0,countCT=0;
			//int revG=0,revGA=0;
			//start to read batman hit file
			int tembuf=1000;
			char Buf[tembuf],Meth[tembuf],Chrom[tembuf], noChrom[tembuf], Strand,context[tembuf],effCT[tembuf];
			printf("\nLoading methratio file...\n");
			FILE* INFILE;//=File_Open(methInfileName,"r");
		for(int f=InFileStart;f<=InFileEnd;f++)
		{
			INFILE=File_Open(argv[f],"r");
			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);
			Init_Progress();
			fgets(Buf,BATBUF,INFILE);//read first header marker..
			while (!feof(INFILE)) 
			{
				Total_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(INFILE);
					off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
					Number_of_Tags=(File_Size/Average_Length)/20;
					Progress=0;
					Show_Progress(Current_Pos*100/File_Size);
				}
				fgets(Meth,BATBUF,INFILE);
				//printf("%s\n", Meth);
//				sscanf(Meth,"%s%u%s%s%*s%*s%d%d",Chrom,&pos,&Strand,context,&countC,&countCT);//,&revG,&revGA
				sscanf(Meth,"%s%u%s%s%d%d",Chrom,&pos,&Strand,context,&countC,&countCT);
				//printf("111\n");	
				if(countCT<Coverage) continue;
				if(countCT > maxcover) continue;
				if(countC > countCT) printf("Wrong pos %d, %d %d\n", pos, countC, countCT);
				int H=String_Hash[Chrom];
				pos--;//0-- for array start 0
				//int whichBins = (int)floor((double)pos/binsStep);//1000
				if(Strand=='+')
				{
					Methy_List[H].plusMethContext[pos]=Context_Hash[context];
					Methy_List[H].plusCount_C[pos]=countC;
					Methy_List[H].plusCount_CT[pos]=countCT;
					//Methy_List[H].binsPlusCount_C[whichBins]+=countC;
					//Methy_List[H].binsPlusCount_CT[whichBins]+=countCT;
				}else if(Strand=='-')
				{
					Methy_List[H].NegMethContext[pos]=Context_Hash[context];
					Methy_List[H].NegCount_C[pos]=countC;
					Methy_List[H].NegCount_CT[pos]=countCT;
					//Methy_List[H].binsNegCount_C[whichBins]+=countC;
					//Methy_List[H].binsNegCount_CT[whichBins]+=countCT;
				}else printf("wrong strand...\n");
				
			}//end read file while
			Done_Progress();
		} //end foreach meth file
			printf("Methratio file Loaded.\n");
//
			Methy_Gff methGff_List[3];
			longerlen = nLevel>bodynLevel?nLevel:bodynLevel;
			printf("%d, %d\n", nLevel, longerlen);
			for ( int i=0;i<3;i++)
			{//0--UP 1--BODY 2--DOWN
				methGff_List[i].CG_C=new long[longerlen]();
				methGff_List[i].CHG_C=new long[longerlen]();
				methGff_List[i].CHH_C=new long[longerlen]();
				methGff_List[i].CG_CT=new long[longerlen]();
				methGff_List[i].CHG_CT=new long[longerlen]();
				methGff_List[i].CHH_CT=new long[longerlen]();
				methGff_List[i].AverPerCG=0;
				methGff_List[i].AverPerCHG=0;
				methGff_List[i].AverPerCHH=0;
				methGff_List[i].AverCG=0;
				methGff_List[i].AverCHG=0;
				methGff_List[i].AverCHH=0;
				methGff_List[i].countCG=0;
				methGff_List[i].countCHG=0;
				methGff_List[i].countCHH=0;
				methGff_List[i].Index=i;
			}
			Methy_Gff methGff_TSS;
			methGff_TSS.CG_C=new long[nLevel]();
			methGff_TSS.CHG_C=new long[nLevel]();
			methGff_TSS.CHH_C=new long[nLevel]();
			methGff_TSS.CG_CT=new long[nLevel]();
			methGff_TSS.CHG_CT=new long[nLevel]();
			methGff_TSS.CHH_CT=new long[nLevel]();

			string du[3];
			du[0]="UP";du[1]="BODY";du[2]="DOWN";
			int FileS=0,FileE=0;
			if(InputGff)
			{
				FileS=GffInFileStart;
				FileE=GffInFileEnd;
			}else if(InputBed)
			{
				FileS=BedInFileStart;
				FileE=BedInFileEnd;
			}else
			{
				fprintf(stderr,"\nNo input Gff/Bed file .\n");
				exit(0);
			}
			for(int f=FileS,ip=0;f<=FileE;f++,ip++){
				AverPerCG=0;AverPerCHG=0;AverPerCHH=0;AverCG=0;AverCHG=0;AverCHH=0;
				printf("\nProcessing %d out of %d. InFile: %s\n", f-FileS+1,FileE-FileS+1, argv[f]);
					FILE* GFFINFILE=File_Open(argv[f],"r");
					
					char T[100];
					//gene body && Promoter
					char G_c[100];
					char G_cg[100];
					char G_chg[100];
					char G_chh[100];
                    char P_c[100];
					char P_cg[100];
					char P_chg[100];
					char P_chh[100];

					char TSS_cg[100];
					char TSS_chg[100];
					char TSS_chh[100];
					//char TSS_cg_matrix[100];
					//char TSS_chg_matrix[100];
					//char TSS_chh_matrix[100];
					char TTS_cg[100];
					char TTS_chg[100];
					char TTS_chh[100];
					char GENE_cg[100];
					char GENE_chg[100];
					char GENE_chh[100];
					char N[100];
					char D[100];
					char TSSprofile[100];
					if(PrefixStart == 0 || ( (PrefixEnd-PrefixStart) != (FileE-FileS) )){
						//only one inputfile or not same
						if(PrefixStart!=0){
							printf("%d %d %d %d\n",PrefixEnd, PrefixStart,FileE, FileS );
							fprintf(stderr, "Warning: The number of prefix is not same as bed/gtf files!! Then the all results will write to one file (%s)!", Output_Name);
						}
						//, f-FileS+1
						sprintf (T, "%s.Methylevel.txt", Output_Name);
						sprintf (G_c, "%s.body.c.txt", Output_Name);
						sprintf (G_cg, "%s.body.cg.txt", Output_Name);
						sprintf (G_chg, "%s.body.chg.txt", Output_Name);
						sprintf (G_chh, "%s.body.chh.txt", Output_Name);
						sprintf (P_c, "%s.Promoter.c.txt", Output_Name);
						sprintf (P_cg, "%s.Promoter.cg.txt", Output_Name);
						sprintf (P_chg, "%s.Promoter.chg.txt", Output_Name);
						sprintf (P_chh, "%s.Promoter.chh.txt", Output_Name);
						//gene
						sprintf (TSS_cg, "%s.TSS.cg.txt", Output_Name);
						sprintf (TSS_chg, "%s.TSS.chg.txt", Output_Name);
						sprintf (TSS_chh, "%s.TSS.chh.txt", Output_Name);
						//sprintf (TSS_cg_matrix, "%s.TSS.cg.matrix", Output_Name);
                	    //sprintf (TSS_chg_matrix, "%s.TSS.chg.matrix", Output_Name);
                	    //sprintf (TSS_chh_matrix, "%s.TSS.chh.matrix", Output_Name);
						sprintf (TTS_cg, "%s.TTS.cg.txt", Output_Name);
						sprintf (TTS_chg, "%s.TTS.chg.txt", Output_Name);
						sprintf (TTS_chh, "%s.TTS.chh.txt", Output_Name);
						//GENE
						sprintf (GENE_cg, "%s.GENE.cg.txt", Output_Name);
						sprintf (GENE_chg, "%s.GENE.chg.txt", Output_Name);
						sprintf (GENE_chh, "%s.GENE.chh.txt", Output_Name);
						//---------------------------------------------------------------------------
						sprintf (N, "%s.AverMethylevel.txt", Output_Name);
						sprintf (D, "%s.annoDensity.txt", Output_Name);
						sprintf (TSSprofile, "%s.TSSprofile.txt", Output_Name);
					}else{
						int Nprefix = PrefixStart+ip;
						sprintf (T, "%s.Methylevel.txt", argv[Nprefix]);
						sprintf (G_c, "%s.body.c.txt", argv[Nprefix]);
						sprintf (G_cg, "%s.body.cg.txt",  argv[Nprefix]);
						sprintf (G_chg, "%s.body.chg.txt", argv[Nprefix]);
						sprintf (G_chh, "%s.body.chh.txt", argv[Nprefix]);
						sprintf (P_c, "%s.Promoter.c.txt", argv[Nprefix]);
						sprintf (P_cg, "%s.Promoter.cg.txt",  argv[Nprefix]);
						sprintf (P_chg, "%s.Promoter.chg.txt", argv[Nprefix]);
						sprintf (P_chh, "%s.Promoter.chh.txt", argv[Nprefix]);
						//gene
						sprintf (TSS_cg, "%s.TSS.cg.txt", argv[Nprefix]);
						sprintf (TSS_chg, "%s.TSS.chg.txt", argv[Nprefix]);
						sprintf (TSS_chh, "%s.TSS.chh.txt", argv[Nprefix]);
						//sprintf (TSS_cg_matrix, "%s.TSS.cg.matrix", argv[Nprefix]);
                	    //sprintf (TSS_chg_matrix, "%s.TSS.chg.matrix", argv[Nprefix]);
                	    //sprintf (TSS_chh_matrix, "%s.TSS.chh.matrix", argv[Nprefix]);
						sprintf (TTS_cg, "%s.TTS.cg.txt", argv[Nprefix]);
						sprintf (TTS_chg, "%s.TTS.chg.txt", argv[Nprefix]);
						sprintf (TTS_chh, "%s.TTS.chh.txt", argv[Nprefix]);
						//GENE
						sprintf (GENE_cg, "%s.GENE.cg.txt", argv[Nprefix]);
						sprintf (GENE_chg, "%s.GENE.chg.txt", argv[Nprefix]);
						sprintf (GENE_chh, "%s.GENE.chh.txt", argv[Nprefix]);
						//---------------------------------------------------------------------------
						sprintf (N, "%s.AverMethylevel.txt", argv[Nprefix]);
						sprintf (D, "%s.annoDensity.txt", argv[Nprefix]);
						sprintf (TSSprofile, "%s.TSSprofile.txt", argv[Nprefix]);
					}

					char mode[10];
					if(PrefixStart == 0 || ( (PrefixEnd-PrefixStart) != (FileE-FileS) ))
						strcpy(mode, "aw+");
					if(f == FileS) 
						strcpy(mode, "w");
					//FILE*
					FILE* MethGff_TSS_OUTFILEcg;
					FILE* MethGff_TSS_OUTFILEchg;
					FILE* MethGff_TSS_OUTFILEchh;
					//FILE* MethGff_TSS_OUTFILEcg_matrix;
                    //FILE* MethGff_TSS_OUTFILEchg_matrix;
                    //FILE* MethGff_TSS_OUTFILEchh_matrix;
					FILE* MethGff_TTS_OUTFILEcg;
					FILE* MethGff_TTS_OUTFILEchg;
					FILE* MethGff_TTS_OUTFILEchh;
					FILE* MethGff_GENE_OUTFILEcg;
					FILE* MethGff_GENE_OUTFILEchg;
					FILE* MethGff_GENE_OUTFILEchh;
					//
					FILE* MethGffbodyOUTFILEc;
					FILE* MethGffbodyOUTFILEcg;
					FILE* MethGffbodyOUTFILEchg;
					FILE* MethGffbodyOUTFILEchh;
					if(Diff)
					{
						MethGffbodyOUTFILEc=File_Open(G_c,mode);
						MethGffbodyOUTFILEcg=File_Open(G_cg,mode);
						MethGffbodyOUTFILEchg=File_Open(G_chg,mode);
						MethGffbodyOUTFILEchh=File_Open(G_chh,mode);
					}
                    FILE* MethGffpromoterOUTFILEc;
					FILE* MethGffpromoterOUTFILEcg;
					FILE* MethGffpromoterOUTFILEchg;
					FILE* MethGffpromoterOUTFILEchh;
					if(PU){
                        MethGffpromoterOUTFILEc=File_Open(P_c,mode);
						MethGffpromoterOUTFILEcg=File_Open(P_cg,mode);
						MethGffpromoterOUTFILEchg=File_Open(P_chg,mode);
						MethGffpromoterOUTFILEchh=File_Open(P_chh,mode);
					}
					//heatmap
					//TSS
					if(TSS)
					{
						MethGff_TSS_OUTFILEcg=File_Open(TSS_cg,mode);
						MethGff_TSS_OUTFILEchg=File_Open(TSS_chg,mode);
						MethGff_TSS_OUTFILEchh=File_Open(TSS_chh,mode);
						//MethGff_TSS_OUTFILEcg_matrix=File_Open(TSS_cg_matrix,mode);
						//MethGff_TSS_OUTFILEchg_matrix=File_Open(TSS_chg_matrix,mode);
						//MethGff_TSS_OUTFILEchh_matrix=File_Open(TSS_chh_matrix,mode);
					}
					//TTS
					if(TTS)
					{
						MethGff_TTS_OUTFILEcg=File_Open(TTS_cg,mode);
						MethGff_TTS_OUTFILEchg=File_Open(TTS_chg,mode);
						MethGff_TTS_OUTFILEchh=File_Open(TTS_chh,mode);
					}
					if(GENE)
					{
						MethGff_GENE_OUTFILEcg=File_Open(GENE_cg,mode);
						MethGff_GENE_OUTFILEchg=File_Open(GENE_chg,mode);
						MethGff_GENE_OUTFILEchh=File_Open(GENE_chh,mode);
					}
					//-----------------------------------------------------------------------------------------------
					printf("\n------------------------------------------------------------------------------------------------\n");
					printf("\nmethylation level output file: %s\n",T);
					printf("Average Methylation level output file: %s\n",N);
					printf("Gff/Bed density output file: %s\n",D);
					if(Diff) printf("Gff/Bed different analysis methylation file: %s.regions.context.xx.txt\n",Output_Name);
					printf("\n------------------------------------------------------------------------------------------------\n");
					//--------------------------------------------------------------------------------------------------
					FILE* MethGffAverOUTFILE=File_Open(T,mode);
				//if(f==FileS) 
				//fprintf(MethGffAverOUTFILE,"#Context\tRegions\tDNA methylation level\n");

				int Progress=0;int Number_of_Tags=INITIAL_PROGRESS_READS;
				fseek(GFFINFILE, 0L, SEEK_END);off64_t File_Size=ftello64(GFFINFILE);rewind(GFFINFILE);//
				Init_Progress();Total_Reads=0;//
				char Gff[BATBUF],temp[BATBUF];
				string id="";
				char Symbol[BATBUF];
				unsigned start=0,end=0;
				Strand = '.';
				float binspan_b = binspan;
				while(fgets(Gff,BATBUF,GFFINFILE)!=0) //(!feof(GFFINFILE))
				{
                    if(Gff[0]=='#') continue;
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
					//fgets(Gff,BATBUF,GFFINFILE);
					//string id="";
					if(InputGff)
					{
						//if(GTF) sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s \"%[^\"]",Chrom,temp,temp,&start,&end,temp,&Strand,temp,temp,Symbol);
			                        //else sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%[^\n\t]",Chrom,temp,temp,&start,&end,temp,&Strand,temp,Symbol);
			            if(GTF) sscanf(Gff,"%s\t%s\t%s\t%d\t%d\t%*s\t%c\t%*s\t%*s \"%[^\"]",Chrom,temp,temp,&start,&end,&Strand,Symbol);
	 			        else    sscanf(Gff,"%s\t%s\t%s\t%u\t%u\t%*s\t%c\t%*s\t%[^;\n\t]",Chrom,temp,temp,&start,&end,&Strand,Symbol);
						for(int i=0;i<strlen(Symbol); i++){
							if(Symbol[i] == ';') {
								Symbol[i]='\0';
								break;
							}
						}
						id = Symbol;
					}else if(InputBed)
					{
						if(bed4){
							sscanf(Gff,"%s\t%u\t%u\t%c",Chrom,&start,&end,&Strand);
							if(showbed == 0) 
							{	
								fprintf(stderr, "bed format and defined strand %c\n", Strand);
								showbed = 1;
							}
						}else if(bed5){
							sscanf(Gff,"%s\t%u\t%u\t%s\t%c",Chrom,&start,&end,Symbol,&Strand);
							if(showbed == 0) 
							{	
								fprintf(stderr, "bed format and defined strand %c\n", Strand);
								showbed = 1;
							}
							id = Symbol;
						}else{
							sscanf(Gff,"%s\t%u\t%u",Chrom,&start,&end);
							Strand='.';
							if(showbed == 0){
								fprintf(stderr, "bed format without strand, if you want define strand, please use bed4 format.\n");
								showbed = 1;
							}
						}
					}
					if(!strcmp(noChrom, Chrom)) continue;
					if(end-start<beddistance) continue;
					if(!InputBed && end-start < beddistance) fprintf(stderr, "Waring: element %s is too shoter to caculate heatmap.\n", id.c_str());
					else if(end-start < beddistance) fprintf(stderr, "Waring: element %s:%d-%d is too shoter to caculate heatmap.\n", Chrom, start, end);
                        	        map<string, int>::iterator it= String_Hash.find(Chrom);
                	                if(it == String_Hash.end()) {
        	                                printf("%s not detected meth\n", Chrom);
	                                        strcpy(noChrom,Chrom);
                                        	continue;
                                	}
					int H=String_Hash[Chrom];
					
					if(start>distance && end+distance <Genome_Offsets[H+1].Offset ) //&& end-start > nLevel 
					{
						start--;end--;
						if(Strand=='-')
						{//Down
							caculate(start-distance,start,Methy_List[H],Strand,methGff_List[2]);
							if(TTS)
							{
								binspan = binspan/2.0;
								caculateHeatmap("TTS",start-distanceHeatmap,start+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_TTS_OUTFILEcg,MethGff_TTS_OUTFILEchg,MethGff_TTS_OUTFILEchh,Chrom,start,end); //, MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
								binspan = binspan_b;
							}
						}else 
						{//Up
							caculate(start-distance,start,Methy_List[H],Strand,methGff_List[0]);
							if(TSS)
							{
								binspan = binspan/2.0;
								caculateHeatmap("TSS",start-distanceHeatmap,start+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_TSS_OUTFILEcg,MethGff_TSS_OUTFILEchg,MethGff_TSS_OUTFILEchh,Chrom,start,end); //, MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
								binspan = binspan_b;
							}
							if(GENE) {
								caculateHeatmap("GENE",start-distanceHeatmap,start,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom,start,end, true); //,MethGff_TSS_OUTFILEcg_matrix, MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
							}
						}
						//body
						
						if(bodyLen != distance){
							binspan = bodyspan;
						}
						caculate(start,end,Methy_List[H],Strand,methGff_List[1],GeneD_List[H]);
						binspan = binspan_b;
						//tss extend 2k，基因内部
						if(Strand=='-')
							caculate(end-distance,end,Methy_List[H],Strand,methGff_TSS);
						else
						{
							caculate(start,start+distance,Methy_List[H],Strand,methGff_TSS);
						}
						
						if(Strand!='-' && GENE){
							//binspan = binspan_b/2;
							if(bodyLen != distance){
								binspan = bodyspan;
							}
							caculateHeatmap("GENE",start,end,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom,start,end, false); //MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix,
							binspan = binspan_b;
						}
						//down
						if(Strand=='-')
						{
							caculate(end,end+distance,Methy_List[H],Strand,methGff_List[0]);
							if(TSS)
							{
								binspan = binspan/2.0;
								caculateHeatmap("TSS",end-distanceHeatmap,end+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_TSS_OUTFILEcg,MethGff_TSS_OUTFILEchg,MethGff_TSS_OUTFILEchh,Chrom,start,end); //, MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
								binspan = binspan_b;
							}
							if(GENE) {
								caculateHeatmap("GENE",end,end+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom, start,end, true); //MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix, 
							}
						}else
						{//down
							caculate(end,end+distance,Methy_List[H],Strand,methGff_List[2]);
							if(TTS)
							{
								binspan = binspan/2.0;
								caculateHeatmap("TTS",end-distanceHeatmap,end+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_TTS_OUTFILEcg,MethGff_TTS_OUTFILEchg,MethGff_TTS_OUTFILEchh,Chrom,start,end); //, MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
								binspan = binspan_b;
							}
							if(GENE) {
								caculateHeatmap("GENE",end, end+distanceHeatmap,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom,start,end, false); //, MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix
							}
						}
						if(Strand=='-' && GENE) {
							//binspan = binspan_b/2;
							if(bodyLen != distance){
								binspan = bodyspan;
							}
							caculateHeatmap("GENE",start,end,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom,start,end,false); // MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix,
							binspan = binspan_b;
							caculateHeatmap("GENE",start-distanceHeatmap,start,Methy_List[H],Strand,id.c_str(),MethGff_GENE_OUTFILEcg,MethGff_GENE_OUTFILEchg,MethGff_GENE_OUTFILEchh,Chrom,start,end, false); //MethGff_TSS_OUTFILEcg_matrix,MethGff_TSS_OUTFILEchg_matrix,MethGff_TSS_OUTFILEchh_matrix,
						}
						if(GENE) {
							fprintf(MethGff_GENE_OUTFILEcg, "\n");
							fprintf(MethGff_GENE_OUTFILEchg, "\n");
							fprintf(MethGff_GENE_OUTFILEchh, "\n");
						}
					}
					//if(Strand=='.') Strand='+';
					if(InputBed && !bed5)
					{
						string temp="";
						start++;end++;
						temp=unsigned2str(start);
						id=Chrom;id+=("."+temp);
						temp=unsigned2str(end);
						id=id+"."+temp;
						start--;end--;
					}
					for(int i=0;i<3;i++)
					{
						if(i==0 && PU)//Promoter
						{
							if(methGff_List[i].AverCG+methGff_List[i].AverCHG+methGff_List[i].AverCHH>0 && methGff_List[i].countCG+methGff_List[i].countCHG+methGff_List[i].countCHH >= nC) fprintf(MethGffpromoterOUTFILEc,"%s\t%d\t%c\tC\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCG+methGff_List[i].AverPerCHG+methGff_List[i].AverPerCHH,methGff_List[i].AverCG+methGff_List[i].AverCHG+methGff_List[i].AverCHH,id.c_str());
							if(methGff_List[i].AverCG>0 && methGff_List[i].countCG >= nC) fprintf(MethGffpromoterOUTFILEcg,"%s\t%d\t%c\tCG\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCG,methGff_List[i].AverCG,id.c_str());
							if(methGff_List[i].AverCHG>0 && methGff_List[i].countCHG >= nC) fprintf(MethGffpromoterOUTFILEchg,"%s\t%d\t%c\tCHG\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCHG,methGff_List[i].AverCHG,id.c_str());
							if(methGff_List[i].AverCHH>0 && methGff_List[i].countCHH >= nC) fprintf(MethGffpromoterOUTFILEchh,"%s\t%d\t%c\tCHH\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCHH,methGff_List[i].AverCHH,id.c_str());
						}
						if(i==1 && Diff)//body
						{
							if(methGff_List[i].AverCG+methGff_List[i].AverCHG+methGff_List[i].AverCHH>0 && methGff_List[i].countCG+methGff_List[i].countCHG+methGff_List[i].countCHH >= nC) fprintf(MethGffbodyOUTFILEc,"%s\t%d\t%c\tC\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCG+methGff_List[i].AverPerCHG+methGff_List[i].AverPerCHH,methGff_List[i].AverCG+methGff_List[i].AverCHG+methGff_List[i].AverCHH,id.c_str());
							if(methGff_List[i].AverCG>0 && methGff_List[i].countCG >= nC) fprintf(MethGffbodyOUTFILEcg,"%s\t%d\t%c\tCG\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCG,methGff_List[i].AverCG,id.c_str());
							if(methGff_List[i].AverCHG>0 && methGff_List[i].countCHG >= nC) fprintf(MethGffbodyOUTFILEchg,"%s\t%d\t%c\tCHG\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCHG,methGff_List[i].AverCHG,id.c_str());
							if(methGff_List[i].AverCHH>0 && methGff_List[i].countCHH >= nC) fprintf(MethGffbodyOUTFILEchh,"%s\t%d\t%c\tCHH\t%d\t%d\t%s\n",Chrom,start+1,Strand,methGff_List[i].AverPerCHH,methGff_List[i].AverCHH,id.c_str());
						}
						if(methGff_List[i].AverCG>0 && methGff_List[i].countCG >= nC) fprintf(MethGffAverOUTFILE,"CG\t%s\t%f\n",du[i].c_str(),(double)methGff_List[i].AverPerCG/(double)methGff_List[i].AverCG);
						if(methGff_List[i].AverCHG>0 && methGff_List[i].countCHG >= nC) fprintf(MethGffAverOUTFILE,"CHG\t%s\t%f\n",du[i].c_str(),(double)methGff_List[i].AverPerCHG/(double)methGff_List[i].AverCHG);
						if(methGff_List[i].AverCHH>0 && methGff_List[i].countCHH >= nC) fprintf(MethGffAverOUTFILE,"CHH\t%s\t%f\n",du[i].c_str(),(double)methGff_List[i].AverPerCHH/(double)methGff_List[i].AverCHH);
					}
					//
				}
				fclose(MethGffAverOUTFILE);
				Done_Progress();

				fprintf(stderr, "\nPrint methresult to files!\n");
				FILE* DensityGffOUTFILE=File_Open(D,mode);
				for ( int i=0;i<Genome_Count;i++)
				{
					unsigned pnCoverage=0,plusCoverage=0,NegCoverage=0,pnCoverage_1=0,plusCoverage_1=0,NegCoverage_1=0;int nbins=0;
					unsigned len=ceil(Genome_Offsets[i+1].Offset/double(chromStep))-1;
					for(int l=0;l<Genome_Offsets[i+1].Offset;l++)
					{
						pnCoverage+=GeneD_List[i].PN_Cover[l];
						plusCoverage+=GeneD_List[i].plusCover[l];
						NegCoverage+=GeneD_List[i].NegCover[l];
				                if( (nbins!=len && l == ((nbins+1)*chromStep)) ||  l==end){
				                    if(nbins<=len &&nbins>0){
				                    	fprintf(DensityGffOUTFILE,"%s\t%d\t%f\t+-\n",Genome_Offsets[i].Genome,nbins-1,(pnCoverage+pnCoverage_1)/(double(chromStep*2)));
								  fprintf(DensityGffOUTFILE,"%s\t%d\t%f\t+\n",Genome_Offsets[i].Genome,nbins-1,(plusCoverage+plusCoverage_1)/(double(chromStep*2)));
								  fprintf(DensityGffOUTFILE,"%s\t%d\t-%f\t-\n",Genome_Offsets[i].Genome,nbins-1,(NegCoverage+NegCoverage_1)/(double(chromStep*2)));
				                    }
				                    plusCoverage_1=plusCoverage,NegCoverage_1=NegCoverage,pnCoverage_1=pnCoverage;
				                    nbins++;
				                    plusCoverage=0,NegCoverage=0,pnCoverage=0;
				                }
					}
				}
				fclose(DensityGffOUTFILE);
				
				FILE* MethGffOUTFILE=File_Open(N,mode);
                                //C
                                //fprintf(MethGffOUTFILE,"C");
                                //for(int i=0;i<3;i++)
                                //{
                                //        for(int j=0;j<nLevel;j++)
                                //        {
                                //                fprintf(MethGffOUTFILE,"\t%f",((double)methGff_List[i].CG_C[j] + (double)methGff_List[i].CHG_C[j] + (double)methGff_List[i].CHH_C[j]) /((double)methGff_List[i].CG_CT[j] + (double)methGff_List[i].CHG_CT[j] + (double)methGff_List[i].CHH_CT[j]));
                                //        }
                                //}
				//CG
				fprintf(MethGffOUTFILE,"CG");
				unsigned printlen = nLevel;
				for(int i=0;i<3;i++)
				{
					if(i==1) printlen = bodynLevel;
					else printlen = nLevel;
					for(int j=0;j<printlen;j++)
					{
						if(methGff_List[i].CG_C[j] > methGff_List[i].CG_CT[j])
							printf("\nError bin %d %d, %ld %ld", i, j, methGff_List[i].CG_C[j], methGff_List[i].CG_CT[j]);
						fprintf(MethGffOUTFILE,"\t%f",(double)methGff_List[i].CG_C[j]/(double)methGff_List[i].CG_CT[j]);
					}
				}
				//CHG
				fprintf(MethGffOUTFILE,"\nCHG");
				for(int i=0;i<3;i++)
				{
					if(i==1) printlen = bodynLevel;
					else printlen = nLevel;
					for(int j=0;j<printlen;j++)
					{
						fprintf(MethGffOUTFILE,"\t%f",(double)methGff_List[i].CHG_C[j]/(double)methGff_List[i].CHG_CT[j]);
					}
				}
				//CHH
				fprintf(MethGffOUTFILE,"\nCHH");
				for(int i=0;i<3;i++)
				{
					if(i==1) printlen = bodynLevel;
					else printlen = nLevel;
					for(int j=0;j<printlen;j++)
					{
						fprintf(MethGffOUTFILE,"\t%f",(double)methGff_List[i].CHH_C[j]/(double)methGff_List[i].CHH_CT[j]);
					}
				}
				fprintf(MethGffOUTFILE,"\n");
				fclose(MethGffOUTFILE);
				//TSS distribution
				FILE* fTSSprofile=File_Open(TSSprofile,mode);
				//CG
				Methy_Gff Meth_out;
				fprintf(fTSSprofile,"CG");
				for(int i=0;i<2;i++)
				{
					if(i==0) Meth_out = methGff_List[i];
					else Meth_out = methGff_TSS;
					for(int j=0;j<nLevel;j++)
					{
						if(Meth_out.CG_C[j] > Meth_out.CG_CT[j])
							printf("\nError bin %d, %ld %ld", i, Meth_out.CG_C[j], Meth_out.CG_CT[j]);
						fprintf(fTSSprofile,"\t%f",(double)Meth_out.CG_C[j]/(double)Meth_out.CG_CT[j]);
					}
				}
				//CHG
				fprintf(fTSSprofile,"\nCHG");
				for(int i=0;i<2;i++)
				{
					if(i==0) Meth_out = methGff_List[i];
					else Meth_out = methGff_TSS;
					for(int j=0;j<nLevel;j++)
					{
						fprintf(fTSSprofile,"\t%f",(double)Meth_out.CHG_C[j]/(double)Meth_out.CHG_CT[j]);
					}
				}
				//CHH
				fprintf(fTSSprofile,"\nCHH");
				for(int i=0;i<2;i++)
				{
					if(i==0) Meth_out = methGff_List[i];
					else Meth_out = methGff_TSS;
					for(int j=0;j<nLevel;j++)
					{
						fprintf(fTSSprofile,"\t%f",(double)Meth_out.CHH_C[j]/(double)Meth_out.CHH_CT[j]);
					}
				}
				fprintf(fTSSprofile,"\n");
				fclose(fTSSprofile);
				
				printf("\n%s\nC : %f\n",argv[f],(double)(AverPerCG+AverPerCHG+AverPerCHH)/(double)(AverCG+AverCHG+AverCHH) );//AverPerCHG
				printf("CG : %f\n",(double)AverPerCG/(double)AverCG);
				printf("CHG : %f\n",(double)AverPerCHG/(double)AverCHG);
				printf("CHH : %f\n",(double)AverPerCHH/(double)AverCHH);
			}
			//delete 
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				delete[] Methy_List[i].plusMethContext;
				delete[] Methy_List[i].NegMethContext;
				delete[] Methy_List[i].plusCount_C ;
				delete[] Methy_List[i].plusCount_CT ;
				delete[] Methy_List[i].NegCount_C;
				delete[] Methy_List[i].NegCount_CT;
				
				delete[] GeneD_List[i].PN_Cover;
				delete[] GeneD_List[i].plusCover;
				delete[] GeneD_List[i].NegCover;
			}
			for ( int i=0;i<3;i++)
			{//0--UP 1--BODY 2--DOWN
				delete[] methGff_List[i].CG_C;
				delete[] methGff_List[i].CHG_C;
				delete[] methGff_List[i].CHH_C;
				delete[] methGff_List[i].CG_CT;
				delete[] methGff_List[i].CHG_CT;
				delete[] methGff_List[i].CHH_CT;
			}
			/*
			delete[] methGff_TSS.CG_C;
			delete[] methGff_TSS.CHG_C;
			delete[] methGff_TSS.CHH_C;
			delete[] methGff_TSS.CG_CT;
			delete[] methGff_TSS.CHG_CT;
			delete[] methGff_TSS.CHH_CT;
			*/
		}
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
		time(&End_Time);printf("\nTime Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	else
		printf("%s \n",Help_String);
}

//{----------------------------------- FILE HANDLING ---------------------------------------------------------

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

//}----------------------------------- FILE HANDLING ---------------------------------------------------------


void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}

void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List,GeneDensity & GeneD_List){
        int step = floor((double(end - start))*binspan);// 5%  2.5% step
		unsigned nLevel=ceil(1/(double)binspan)-1;
        int nbins=0;
        methGff_List.AverPerCG=0;methGff_List.AverCG=0;
        methGff_List.AverPerCHG=0;methGff_List.AverCHG=0;
        methGff_List.AverPerCHH=0;methGff_List.AverCHH=0;
        
            long countperCG=0,countperCHG=0,countperCHH=0,countCG=0,countCHG=0,countCHH=0;
            long countperCG_1=0,countperCHG_1=0,countperCHH_1=0,countCG_1=0,countCHG_1=0,countCHH_1=0;    
            for(int i=start;i<=end;i++)
            {
            	GeneD_List.PN_Cover[i]=1;
				//if(MethyList.plusCount_CT[i] < MethyList.plusCount_C[i]) continue;
	                //context
	                if(Strand=='+' || Strand=='.')
	                {
	                	GeneD_List.plusCover[i]=1;
	                        if( MethyList.plusMethContext[i] ==1 ){
	                            countperCG+= MethyList.plusCount_C[i];
	                            countCG+=MethyList.plusCount_CT[i]; 
	                            methGff_List.AverPerCG +=MethyList.plusCount_C[i];
	                            methGff_List.AverCG +=MethyList.plusCount_CT[i]; 
	                            AverPerCG +=MethyList.plusCount_C[i];
	                            AverCG +=MethyList.plusCount_CT[i]; 
	                            methGff_List.countCG++;
	                        }else if( MethyList.plusMethContext[i] ==2 ){
	                            countperCHG+=MethyList.plusCount_C[i];
	                            countCHG+=MethyList.plusCount_CT[i];
	                            methGff_List.AverPerCHG +=MethyList.plusCount_C[i];
	                            methGff_List.AverCHG += MethyList.plusCount_CT[i];
	                            AverPerCHG +=MethyList.plusCount_C[i];
	                            AverCHG +=MethyList.plusCount_CT[i]; 
	                            methGff_List.countCHG++;
	                        }else if( MethyList.plusMethContext[i] ==3 ){
	                            countperCHH+=MethyList.plusCount_C[i];
	                            countCHH+=MethyList.plusCount_CT[i];
	                            methGff_List.AverPerCHH +=MethyList.plusCount_C[i];
	                            methGff_List.AverCHH +=MethyList.plusCount_CT[i];
	                            AverPerCHH +=MethyList.plusCount_C[i];
	                            AverCHH +=MethyList.plusCount_CT[i]; 
	                            methGff_List.countCHH++;
	                        }
			}
			if(Strand=='-' || Strand=='.')
			{
				GeneD_List.NegCover[i]=1;
	                        if( MethyList.NegMethContext[i] ==1 ){
	                            countperCG+= MethyList.NegCount_C[i];
	                            countCG+=MethyList.NegCount_CT[i]; 
	                            methGff_List.AverPerCG +=MethyList.NegCount_C[i];
	                            methGff_List.AverCG +=MethyList.NegCount_CT[i]; 
	                            AverPerCG +=MethyList.NegCount_C[i];
	                            AverCG +=MethyList.NegCount_CT[i]; 
	                            methGff_List.countCG++;
	                        }else if( MethyList.NegMethContext[i] ==2 ){
	                            countperCHG+=MethyList.NegCount_C[i];
	                            countCHG+=MethyList.NegCount_CT[i];
	                            methGff_List.AverPerCHG +=MethyList.NegCount_C[i];
	                            methGff_List.AverCHG += MethyList.NegCount_CT[i];
	                            AverPerCHG +=MethyList.NegCount_C[i];
	                            AverCHG +=MethyList.NegCount_CT[i]; 
	                            methGff_List.countCHG++;
	                        }else if( MethyList.NegMethContext[i] ==3 ){
	                            countperCHH+=MethyList.NegCount_C[i];
	                            countCHH+=MethyList.NegCount_CT[i];
	                            methGff_List.AverPerCHH +=MethyList.NegCount_C[i];
	                            methGff_List.AverCHH +=MethyList.NegCount_CT[i];
	                            AverPerCHH +=MethyList.NegCount_C[i];
	                            AverCHH +=MethyList.NegCount_CT[i]; 
	                            methGff_List.countCHH++;
	                        }
			}
                if( ( (nbins!=nLevel && (i-start) == ((nbins+1)*step)) ||  i==end ) && end-start > nLevel ){
                	if(i==end) nbins=nLevel;
                    if(nbins<=nLevel &&nbins>0){
                    	if(Strand=='+' || Strand=='.')
                    	{
	                            methGff_List.CG_C[nbins-1] += (countperCG+countperCG_1);
	                            methGff_List.CG_CT[nbins-1] += (countCG+countCG_1);
	                            methGff_List.CHG_C[nbins-1] += (countperCHG+countperCHG_1);
	                            methGff_List.CHG_CT[nbins-1] += (countCHG+countCHG_1);
	                            methGff_List.CHH_C[nbins-1] += (countperCHH+countperCHH_1);
	                            methGff_List.CHH_CT[nbins-1] += (countCHH+countCHH_1);
                        }
                        else
                        {
	                            methGff_List.CG_C[nLevel-nbins] += (countperCG+countperCG_1);
	                            methGff_List.CG_CT[nLevel-nbins] += (countCG+countCG_1);
	                            methGff_List.CHG_C[nLevel-nbins] += (countperCHG+countperCHG_1);
	                            methGff_List.CHG_CT[nLevel-nbins] += (countCHG+countCHG_1);
	                            methGff_List.CHH_C[nLevel-nbins] += (countperCHH+countperCHH_1);
	                            methGff_List.CHH_CT[nLevel-nbins] += (countCHH+countCHH_1);
                        }
						
                    }
                    countperCG_1=countperCG;countperCHG_1=countperCHG;countperCHH_1=countperCHH;countCG_1=countCG;countCHG_1=countCHG;countCHH_1=countCHH;
                    nbins++;
                    countperCG=0;countperCHG=0;countperCHH=0;countCG=0;countCHG=0;countCHH=0;
                    
                }
      }
}
void caculate(int start,int end,Methy_Hash MethyList,char Strand,Methy_Gff & methGff_List){
        int step = floor((double(end - start))*binspan);// 5%  2.5% step
		unsigned nLevel=ceil(1/(double)binspan)-1;
        int nbins=0;
        methGff_List.AverPerCG=0;methGff_List.AverCG=0;
        methGff_List.AverPerCHG=0;methGff_List.AverCHG=0;
        methGff_List.AverPerCHH=0;methGff_List.AverCHH=0;
        
            long countperCG=0,countperCHG=0,countperCHH=0,countCG=0,countCHG=0,countCHH=0;
            long countperCG_1=0,countperCHG_1=0,countperCHH_1=0,countCG_1=0,countCHG_1=0,countCHH_1=0;    
            for(int i=start;i<=end;i++)
            {
				//if(MethyList.plusCount_CT[i] < MethyList.plusCount_C[i])
				//	continue;
	                //context
	               if(Strand=='+' || Strand=='.')
	                {
	                        if( MethyList.plusMethContext[i] ==1 ){
	                            countperCG+= MethyList.plusCount_C[i];
	                            countCG+=MethyList.plusCount_CT[i]; 
	                            methGff_List.AverPerCG +=MethyList.plusCount_C[i];
	                            methGff_List.AverCG +=MethyList.plusCount_CT[i]; 
	                            methGff_List.countCG++;
	                        }else if( MethyList.plusMethContext[i] ==2 ){
	                            countperCHG+=MethyList.plusCount_C[i];
	                            countCHG+=MethyList.plusCount_CT[i];
	                            methGff_List.AverPerCHG +=MethyList.plusCount_C[i];
	                            methGff_List.AverCHG += MethyList.plusCount_CT[i];
	                            methGff_List.countCHG++;
	                        }else if( MethyList.plusMethContext[i] ==3 ){
	                            countperCHH+=MethyList.plusCount_C[i];
	                            countCHH+=MethyList.plusCount_CT[i];
	                            methGff_List.AverPerCHH +=MethyList.plusCount_C[i];
	                            methGff_List.AverCHH +=MethyList.plusCount_CT[i];
	                            methGff_List.countCHH++;
	                        }
			}
			if(Strand=='-' || Strand=='.')
			{
	                        if( MethyList.NegMethContext[i] ==1 ){
	                            countperCG+= MethyList.NegCount_C[i];
	                            countCG+=MethyList.NegCount_CT[i]; 
	                            methGff_List.AverPerCG +=MethyList.NegCount_C[i];
	                            methGff_List.AverCG +=MethyList.NegCount_CT[i]; 
	                            methGff_List.countCG++;
	                        }else if( MethyList.NegMethContext[i] ==2 ){
	                            countperCHG+=MethyList.NegCount_C[i];
	                            countCHG+=MethyList.NegCount_CT[i];
	                            methGff_List.AverPerCHG +=MethyList.NegCount_C[i];
	                            methGff_List.AverCHG += MethyList.NegCount_CT[i];
	                            methGff_List.countCHG++;
	                        }else if( MethyList.NegMethContext[i] ==3 ){
	                            countperCHH+=MethyList.NegCount_C[i];
	                            countCHH+=MethyList.NegCount_CT[i];
	                            methGff_List.AverPerCHH +=MethyList.NegCount_C[i];
	                            methGff_List.AverCHH +=MethyList.NegCount_CT[i];
	                            methGff_List.countCHH++;
	                        }
			}
                if( ((nbins!=nLevel && (i-start) == ((nbins+1)*step)) ||  i==end) && end-start > nLevel){
                	if(i==end) nbins=nLevel;
                    if(nbins<=nLevel &&nbins>0){
                    	if(Strand=='+' || Strand=='.')
                    	{
	                            methGff_List.CG_C[nbins-1] += (countperCG+countperCG_1);
	                            methGff_List.CG_CT[nbins-1] += (countCG+countCG_1);
	                            methGff_List.CHG_C[nbins-1] += (countperCHG+countperCHG_1);
	                            methGff_List.CHG_CT[nbins-1] += (countCHG+countCHG_1);
	                            methGff_List.CHH_C[nbins-1] += (countperCHH+countperCHH_1);
	                            methGff_List.CHH_CT[nbins-1] += (countCHH+countCHH_1);
                        }
                        else 
                        {
	                            methGff_List.CG_C[nLevel-nbins] += (countperCG+countperCG_1);
	                            methGff_List.CG_CT[nLevel-nbins] += (countCG+countCG_1);
	                            methGff_List.CHG_C[nLevel-nbins] += (countperCHG+countperCHG_1);
	                            methGff_List.CHG_CT[nLevel-nbins] += (countCHG+countCHG_1);
	                            methGff_List.CHH_C[nLevel-nbins] += (countperCHH+countperCHH_1);
	                            methGff_List.CHH_CT[nLevel-nbins] += (countCHH+countCHH_1);
                        }
                    }
                    countperCG_1=countperCG;countperCHG_1=countperCHG;countperCHH_1=countperCHH;countCG_1=countCG;countCHG_1=countCHG;countCHH_1=countCHH;
                    nbins++;
                    countperCG=0;countperCHG=0;countperCHH=0;countCG=0;countCHG=0;countCHH=0;
                    
                }
      }
}

//caculate gene heatmap //http://www.sciencedirect.com/science/article/pii/S0092867413002225  //Fig.6G
void caculateHeatmap(const char* type,int start,int end,Methy_Hash MethyList,char Strand,const char* id,FILE *methGffcg,FILE *methGffchg,FILE *methGffchh,char* chrom, int geneS, int geneE , bool printtitle){ //FILE *methGffcg_matrix,FILE *methGffchg_matrix,FILE *methGffchh_matrix ,
        int step = floor((double(end - start))*((double)binspan/2)); //5%  2.5% step
	unsigned nLevel=ceil(1/((double)binspan/2))-1;
//	int step = floor((double(end - start))*((double)binspan/1));
//	unsigned nLevel=ceil(1/((double)binspan/1))-1;
        int nbins=0;
        double* Smeth_cg=new double[nLevel]();double* Smeth_chg=new double[nLevel]();double* Smeth_chh=new double[nLevel]();
        
            int countperCG=0,countperCHG=0,countperCHH=0,countCG=0,countCHG=0,countCHH=0;
            int countperCG_1=0,countperCHG_1=0,countperCHH_1=0,countCG_1=0,NegcountCG_1=0,countCHG_1=0,countCHH_1=0;

            for(int i=start;i<=end;i++)
            {
	                //context
	               if(Strand=='+' || Strand=='.')
	                {
	                        if( MethyList.plusMethContext[i] ==1 ){
	                            countperCG+= MethyList.plusCount_C[i];
	                            countCG+=MethyList.plusCount_CT[i];
	                        }else if( MethyList.plusMethContext[i] ==2 ){
	                            countperCHG+=MethyList.plusCount_C[i];
	                            countCHG+=MethyList.plusCount_CT[i];
	                        }else if( MethyList.plusMethContext[i] ==3 ){
	                            countperCHH+=MethyList.plusCount_C[i];
	                            countCHH+=MethyList.plusCount_CT[i];
	                        }
			}
			else if(Strand=='-' )//|| Strand=='.'
			{
	                        if( MethyList.NegMethContext[i] ==1 ){
	                            countperCG+= MethyList.NegCount_C[i];
	                            countCG+=MethyList.NegCount_CT[i]; 
	                        }else if( MethyList.NegMethContext[i] ==2 ){
	                            countperCHG+=MethyList.NegCount_C[i];
	                            countCHG+=MethyList.NegCount_CT[i];
	                        }else if( MethyList.NegMethContext[i] ==3 ){
	                            countperCHH+=MethyList.NegCount_C[i];
	                            countCHH+=MethyList.NegCount_CT[i];
	                        }
			}
                if( ((nbins!=nLevel && (i-start) == ((nbins+1)*step)) ||  i==end) && end-start > nLevel){
                    if(nbins<=nLevel &&nbins>0){
                    	if(Strand=='+' || Strand=='.')
                    	{
					  //if(countCG+countCG_1>0) fprintf(methGffcg,"\t%f",double(countperCG+countperCG_1)/double(countCG+countCG_1));
					  if(countCG+countCG_1>0) Smeth_cg[nbins-1]=double(countperCG+countperCG_1)/(double(countCG+countCG_1));
					  if(countCG+countCG_1<=5) Smeth_cg[nbins-1]=-1;
					  if(countCHG+countCHG_1>0) Smeth_chg[nbins-1]=double(countperCHG+countperCHG_1)/(double(countCHG+countCHG_1));
					  if(countCHG+countCHG_1<=5) Smeth_chg[nbins-1] = -1;
					  if(countCHH+countCHH_1>0) Smeth_chh[nbins-1]=double(countperCHH+countperCHH_1)/(double(countCHH+countCHH_1));
					  if(countCHH+countCHH_1<=5) Smeth_chh[nbins-1]=-1;
                        }
                        else
                        {
					  if(countCG+countCG_1>0) Smeth_cg[nLevel-nbins]=double(countperCG+countperCG_1)/(double(countCG+countCG_1));
					  if(countCG+countCG_1<=5) Smeth_cg[nLevel-nbins] = -1;
					  if(countCHG+countCHG_1>0) Smeth_chg[nLevel-nbins]=double(countperCHG+countperCHG_1)/(double(countCHG+countCHG_1));
					  if(countCHG+countCHG_1<=5) Smeth_chg[nLevel-nbins] = -1;
					  if(countCHH+countCHH_1>0) Smeth_chh[nLevel-nbins]=double(countperCHH+countperCHH_1)/(double(countCHH+countCHH_1));
					  if(countCHH+countCHH_1<=5) Smeth_chh[nLevel-nbins] = -1;
                        }
                    }
                    countperCG_1=countperCG;countperCHG_1=countperCHG;countperCHH_1=countperCHH;countCG_1=countCG;countCHG_1=countCHG;countCHH_1=countCHH;
                    nbins++;
                    countperCG=0;countperCHG=0;countperCHH=0;countCG=0;countCHG=0;countCHH=0;
                    if(i==end) nbins=nLevel;
                }
      }
    if(printtitle && (!strcmp(type,"GENE") || !strcmp(type,"TSS")))
	{
		if(!strcmp(type,"GENE")){
			/*
			fprintf(methGffcg,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			fprintf(methGffchg,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			fprintf(methGffchh,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			*/
			fprintf(methGffcg,"%s",id);
			fprintf(methGffchg,"%s",id);
			fprintf(methGffchh,"%s",id);
		}else{
			//fprintf(methGffcg_matrix,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			//fprintf(methGffchg_matrix,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			//fprintf(methGffchh_matrix,"%s\t%d\t%d\t%s\t.\t.",chrom,geneS+1,geneE+1,id);
			fprintf(methGffcg,"%s",id);
			fprintf(methGffchg,"%s",id);
			fprintf(methGffchh,"%s",id);
		}
	}else if(printtitle) {
		fprintf(methGffcg,"%s",id);
	    fprintf(methGffchg,"%s",id);
		fprintf(methGffchh,"%s",id);
	}
        unsigned cut = ceil((double)nLevel/2);
        unsigned BeginNo_cg=0,BeginNo_chg=0,BeginNo_chh=0;
        //if(!strcmp(type,"TSS"))
        {
	        for(int i=0;i<nLevel;i++)
	        {
	        	if(BeginNo_cg==0 && !strcmp(type,"TSS"))
	        	{	
		        	if(i>=cut )
		       		{
		       		if(Smeth_cg[i]>0) BeginNo_cg=i;
		       		else if(i==nLevel-1) BeginNo_cg=0;
		        	}
	        	}  
	        	if(Smeth_cg[i]>1) Smeth_cg[i]=0;
			if(Smeth_cg[i] < 0){
				if(!strcmp(type, "GENE")) fprintf(methGffcg,"\tnan");
				else if(!strcmp(type, "TSS")) {
					fprintf(methGffcg,"\t0");
					//fprintf(methGffcg_matrix,"\tnan");
				}else fprintf(methGffcg,"\t0");
			}else {
				fprintf(methGffcg,"\t%f",Smeth_cg[i]);
				//fprintf(methGffcg_matrix,"\t%f", Smeth_cg[i]);
			}
	        	if(BeginNo_chg==0 && !strcmp(type,"TSS"))
	        	{	
		        	if(i>=cut)
		       		{
		       		if(Smeth_chg[i]>0) BeginNo_chg=i;
		       		else if(i==nLevel-1) BeginNo_chg=0;
		        	}
	        	}  
	        	if(Smeth_chg[i]>1) Smeth_chg[i]=0;
				if(Smeth_chg[i] < 0){
					if(!strcmp(type, "GENE")) fprintf(methGffchg,"\tnan");
					else if(!strcmp(type, "TSS")) {
						fprintf(methGffchg,"\t0");
						//fprintf(methGffchg_matrix,"\tnan");
					}else fprintf(methGffchg,"\t0");
				}else {
				fprintf(methGffchg,"\t%f",Smeth_chg[i]);
				//fprintf(methGffchg_matrix,"\t%f",Smeth_chg[i]);
			}
	        	if(BeginNo_chh==0 && !strcmp(type,"TSS"))
	        	{	
		        	if(i>=cut)
		       		{
		       		if(Smeth_chh[i]>0) BeginNo_chh=i;
		       		else if(i==nLevel-1) BeginNo_chh=0;
		        	}
	        	}  
	        	if(Smeth_chh[i]>1) Smeth_chh[i]=0;
			if(Smeth_chh[i] < 0){
                                if(!strcmp(type, "GENE")) fprintf(methGffchh,"\tnan");
                                else if(!strcmp(type, "TSS")) {
                                        fprintf(methGffchh,"\t0");
                                        //fprintf(methGffchh_matrix,"\tnan");
                                }else fprintf(methGffchh,"\t0");
                        }else {
				fprintf(methGffchh,"\t%f",Smeth_chh[i]);
				//fprintf(methGffchh_matrix,"\t%f", Smeth_chh[i]);
			}
	        }
        }
        
        if(!strcmp(type,"TTS"))
        {
	        for(int i=nLevel-1;i>=0;i--)
	        {
	        	if(BeginNo_cg==0)
	        	{	
		        	if(i< cut-1 )
		       	{
		       		if(Smeth_cg[i]>0) BeginNo_cg=i;
		       		else if(i==0) BeginNo_cg=cut;
		        	}
	        	}  
	        	if(BeginNo_chg==0)
	        	{	
		        	if(i<cut-1)
		       	{
		       		if(Smeth_chg[i]>0) BeginNo_chg=i;
		       		else if(i==0) BeginNo_chg=cut;
		        	}
	        	}  
	        	if(BeginNo_chh==0)
	        	{	
		        	if(i<cut-1)
		       	{
		       		if(Smeth_chh[i]>0) BeginNo_chh=i;
		       		else if(i==0) BeginNo_chh=cut;
		        	}
	        	}  
	        }
        }
	if(strcmp(type,"GENE")!=0){
        	//fprintf(methGffcg,"\t%d\n",BeginNo_cg);
			fprintf(methGffcg,"\n");
		if(strcmp(type,"TSS")==0) {
			//fprintf(methGffcg_matrix,"\n");
			//fprintf(methGffchg_matrix,"\n");
			//fprintf(methGffchh_matrix,"\n");
		}
	        //fprintf(methGffchg,"\t%d\n",BeginNo_chg);
        	//fprintf(methGffchh,"\t%d\n",BeginNo_chh);
			fprintf(methGffchg,"\n");
        	fprintf(methGffchh,"\n");
	}
}
