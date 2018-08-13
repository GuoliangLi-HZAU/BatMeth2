#include <iostream>
#include <string>
#include <stdio.h>
#include <deque>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>

int chromname_len=400;
const int BUFSIZE=600;
FILE* File_Open(const char* File_Name,const char* Mode);
void rrbs_index(FILE* OUTrrbsfile, char* chr_id, std::string seq, int low_bound, int up_bound, char* redss);
int NLen=50;
inline char ConvToUpperA(char chConv)
{

    return (chConv >= 'a' && chConv <= 'z')? (chConv & 0xdf) : chConv;

}
int main(int argc, char* argv[]){
        fprintf(stderr, "\nBatMeth2: build index rrbs");
        const char* Help_String="\nCommand Format : BatMeth2 build_index [options] genome.fa\n"
                "\nUsage:\n"
                "\t--lowD     lower distance of fragment length, default:20\n"
                "\t--upD      upper distance of fragment length, default:600\n"
                "\t-S         Set restriction enzyme digestion sites. default: [-S C-CGG] for MspI digestion.\n"
                "\t-h|--help";


	//restriction enzyme digestion sites
	char* redss = new char[100];
	strcpy(redss, "C-CGG");
	char* genome_file = new char[100];
	int low_distance = 20;
	int up_distance = 500;


        if(argc<2)
        {
                printf("%s \n",Help_String);
                exit(0);
        }
        for(int i=1;i<argc;i++)
        {
                if(!strcmp(argv[i], "--lowD") )
                {
                        low_distance=atoi(argv[++i]);
                }else if(!strcmp(argv[i], "--upD")){
			up_distance=atoi(argv[++i]);
		}else if(!strcmp(argv[i], "-S")){
			strcpy(redss, argv[++i]);
		}
	}
	strcpy(genome_file, argv[argc-1]);
	
	//PROCESS redss
	int i=0;bool movp=false;
	for(i=0;i<strlen(redss);i++){
		if(redss[i]=='-'){
			if(i+1<strlen(redss))
				redss[i]=redss[i+1];
			else redss[i]=0;
			movp=true;
		}else if(movp){
			if(i+1<strlen(redss))
                        	redss[i]=redss[i+1];
			else redss[i]=0;
		}
	}
	
	fprintf(stderr, "\nrestriction enzyme digestion sites: %s\n", redss);

	std::string seq="";
	//read genome file and process with reds
	FILE* INFILE=File_Open(genome_file,"r");

	char* Buffer=new char[BUFSIZE];
	Buffer[BUFSIZE]=0;
	char* chr_id=new char[chromname_len];
	chr_id[chromname_len]=0;
	char* genome_rrbs=new char[strlen(genome_file)+20];
	strcpy(genome_rrbs, genome_file);
	strcat(genome_rrbs, ".rrbs.fa");

	FILE * OUTrrbs=File_Open(genome_rrbs, "w");

	while (fgets(Buffer,BUFSIZE,INFILE))
	{
		if(Buffer[strlen(Buffer)-1]=='\n')
			Buffer[strlen(Buffer)-1]='\0';
		if (Buffer[0] != '>') 
		{
			seq+=Buffer;
		}
		else if (Buffer[0] == '>') 
		{
			if(seq!=""){
				fprintf(stderr, "\nProcess chromosome: %s", chr_id);
				rrbs_index(OUTrrbs, chr_id, seq, low_distance, up_distance, redss);
			}
			strcpy(chr_id, Buffer);
			seq="";
		}
		
	}
	//the last chrom
	fprintf(stderr, "\nProcess chromosome: %s", chr_id);
	rrbs_index(OUTrrbs, chr_id, seq, low_distance, up_distance, redss);
	fprintf(stderr, "\nDone!\n");
	fclose(INFILE);
	fclose(OUTrrbs);


	delete []Buffer;
	delete []chr_id;
	delete []genome_file;
	delete []genome_rrbs;
	delete []redss;
}

FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		fprintf(stderr, "File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

void rrbs_index(FILE* OUTrrbsfile, char* chr_id, std::string seq, int low_bound, int up_bound, char* redss){
	//transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	unsigned int cLen = seq.length();
        for (size_t i = 0; i<cLen; ++i)
        {
	        seq[i] = ConvToUpperA(seq[i]);
        }

	int redlen=strlen(redss);
    //-- get all red sites ---------------------------------
	fprintf(stderr, "\nDectet restriction enzyme digestion sites in the chromosome.");

	std::deque<std::size_t> reds_loci;
	reds_loci.push_back(0);
	int num=0;
	std::size_t found = seq.find(redss, 0);
	while (found!=std::string::npos){
		num++;
		reds_loci.push_back(found);
		found = seq.find(redss, found+1);
	}
	reds_loci.push_back(cLen-redlen);

	//std::string map_seq=""; //too slow for map_seq = map_seq + 'N'
	char* map_seq = new char[cLen+1];
	strcpy(map_seq, "None");

	fprintf(stderr, "\nProcess the genome with restriction enzyme digestion sites number: %d, genomeLen: %d", reds_loci.size(), cLen);
	unsigned m=0;
	int rstart=0;int rend=0;int distance=0;
	if(reds_loci.size()>1){
		rstart=reds_loci[0];
		rend=reds_loci[1]+redlen;
		distance=rend-rstart;
	}
	bool valid=false;int index=0;
	while (m<cLen){
		if(reds_loci.size()>1){
			if (m >= rstart &&  m < rend-redlen){
				if(distance>=low_bound && distance <= up_bound){
					map_seq[index] = seq[m];
					index++;
					valid=true;
				}else{
					map_seq[index] = '-'; //'N'
					index++;
					valid=false;
				}
			}else{
				if(m<rstart){
					m++;
					continue;
				}
				if(m < rend && valid){
					map_seq[index] = seq[m];
					index++;
					m++; continue;
				}else m--;
				reds_loci.pop_front();
				if(reds_loci.size()>1){
					rstart=reds_loci[0];
					rend=reds_loci[1]+redlen;
					distance=rend-rstart;
				}
			}
		}
		m++;
	}
	map_seq[cLen]='\0';
	if(!strcmp(map_seq, "None")) strcpy(map_seq, seq.c_str());
	fprintf(stderr, "\nModify output length: %s", chr_id);
	char* _map_seq= new char[cLen+1+(int)(cLen/NLen)];
	int i;int nLen=1;int step=0;
	for(i=0;i<cLen;i++,nLen++){
		if(nLen%NLen == 0){
			_map_seq[i+step]=map_seq[i];
			step=(int)(nLen/NLen);
			_map_seq[i+step]='\n';
		}else _map_seq[i+step]=map_seq[i];
	}
	fprintf(stderr, "\nPrint the modified genome.");
	fprintf(OUTrrbsfile, "%s\n%s\n", chr_id, _map_seq);

	delete []map_seq;
	delete []_map_seq;
}

