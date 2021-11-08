#include <iostream>
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <limits.h>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <unistd.h>
#include "zlib.h"
#include <ctype.h> //isdigit
#define MAXDES 500 
#define MAXTAG 1000
#define FQ	2
#define FA	3
using std::string;
using std::vector;
int MAX_PATH = 1000;

//g++ memalign.cpp -o memalign -lz -lpthread
//g++ -fPIE memalign.cpp -o memalign -I./samtools-0.1.18/ -L./samtools-0.1.18/ -lbam -lz -lpthread

#ifdef __cplusplus
extern "C"
{
#endif
        #include "./samtools-0.1.18/sam_header.h"
        #include "./samtools-0.1.18/sam.h"
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif
        bam_header_t *bam_header_dup(const bam_header_t *h0);
        samfile_t *samopen(const char *fn, const char *mode, const void *aux);
        void samclose(samfile_t *fp);
#ifdef __cplusplus
}
#endif

bool PAIRED = true;
bool gzinfile = false;
pthread_mutex_t Lock_gzfile;
bool ALLMAP = false;
//alignment
int threads = 8;
float mismatches = 0.1;
string genome_index = "";
string outformat = "bam";

// fastp
string fastp="";

//aligner
string aligner = "bwame"; // bwa-meth bsmap bismark
string genome_others = "";

string abspath = "";
string workdir;
bool cleanreads=true;
string programname;

// fq format
int gzfilelist[1024];
int gzfilelist_1[1024];
int gzfilelist_2[1024];

void usage(){
    fprintf(stderr, "\nRMAT [mode] [paramaters]\n");
    fprintf(stderr, "mode:  build_index, pipel, align, calmeth, annoation, methyPlot, batDMR, visul2sample, DMCplot\n\n");
    fprintf(stderr, "Example:\n   RMat pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -1 in1.fq.gz -2 in2.fq.gz --gff gene.gff -o outprefix\n");
    fprintf(stderr, "Or single-end:\n   RMat pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -i in.fq.gz --gff gene.gff -o outprefix\n\n");
//    ~/software_devp/batmeth2/bin/test pipel --aligner bwa-meth --go ~/practice/Genome/arabidopsis/arabidopsis_bwa-meth/TAIR10_chr_all.fa --fastp ~/software/fastp/fastp -g ~/practice/Genome/arabidopsis/arabidopsis_batmeth2_index/TAIR10_chr_all.fa -1 test1.fq -2 test2.fq --gff ~/practice/Genome/arabidopsis/TAIR10.gene.modify.gff -f 1 -o pipel.clean
    fprintf(stderr, "\n[build_index]\n");
    fprintf(stderr, "    Usage: rmat build_index genomefile. (wgbs data, must run this step first), or\n");

    fprintf(stderr, "\n[pipel (Contains: align, calmeth, annoation, methyPlot)]\n");
    fprintf(stderr, "[fastp location]\n");
    fprintf(stderr, "    --fastp    fastp program location.\n");
    fprintf(stderr, "               ** If --fastp is not defined, the input file should be clean data.\n");
    fprintf(stderr, "\n[main paramaters]\n");
    fprintf(stderr, "    --config [config file].   When we run pipel function in batches datasets, \n");
    fprintf(stderr, "                              please fill in the specified configuration file. \n");
    fprintf(stderr, "                              And there is a sample file (multirun.onf) in the RMat directory.\n");
    fprintf(stderr, "    --mp [4]                  When batch processing data, we set the number of samples to run at a time (-mp, default is 4), and each sample needs six threads (- P parameter) by default.\n");
    fprintf(stderr, "    --ap [24]                 When batch processing data, the number of samples to run at a time (default 4) * threads of each sample (default 6)\n");
    fprintf(stderr, "    -o [outprefix]            Name of output file prefix\n");
    fprintf(stderr, "    -O [out folder]           Output of result file to specified folder, default output to current folder (./)\n");
    fprintf(stderr, "    -of [SAM/BAM]             Output format, default BAM.\n");
    fprintf(stderr, "\n[alignment paramaters]\n");
    fprintf(stderr, "    -i    Name of input file, support .fq/.fastq and .gz/.gzip format. if paired-end. please use -1, -2\n");
    fprintf(stderr, "    -1    Name of input file left end, if single-end. please use -i\n");
    fprintf(stderr, "    -2    Name of input file right end\n");
    fprintf(stderr, "          -i/-1/-2 can be comma-separated lists (no whitespace), only supported in RMat aligner.\n");
    fprintf(stderr, "    -g    Name of the genome mapped against, make sure build index first.\n");
    fprintf(stderr, "    -p    Launch <integer> threads\n");
    fprintf(stderr, "\n[calmeth paramaters]\n");
    fprintf(stderr, "    --Qual      calculate the methratio while read QulityScore >= Q. default:20\n");
    fprintf(stderr, "    --redup     REMOVE_DUP, 0 or 1, default 1\n");
    fprintf(stderr, "    --region    Bins for DMR calculate , default 1000bp .\n");
    fprintf(stderr, "    -f          for sam format outfile contain methState. [0 or 1], default: 0 (dont output this file).\n");
	fprintf(stderr, "    -n          maximum mismatches allowed due to seq. default 0.1 percentage of the read length. [0-0.3]\n");
    fprintf(stderr, "\n[calmeth and annoation paramaters]\n");
    fprintf(stderr, "    --coverage    >= <INT> coverage. default:4\n");
    fprintf(stderr, "    --binCover    >= <INT> nCs per region. default:1\n");
    fprintf(stderr, "    --chromstep   Chromsome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp)\n");
    fprintf(stderr, "\n[annoation paramaters]\n");
    fprintf(stderr, "    --gtf/--gff/--bed    Gtf or gff file / bed file\n");
    fprintf(stderr, "    --distance           RNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000\n");
    fprintf(stderr, "    --step               Gene body and their flanking sequences using an overlapping sliding window of 0.02 of the sequence length at a step of 0.01 of the sequence length. So default step: 0.01 (1%)\n");
    fprintf(stderr, "    -C                   <= <INT> coverage. default:1000\n");
    fprintf(stderr, "\n[MethyPlot paramaters]\n");
    fprintf(stderr, "    --CG       CG ratio for heatmap, [0-1], default 0.6\n");
    fprintf(stderr, "    --CHG      CHG ratio for heatmap, [0-1], default 0.2\n");
    fprintf(stderr, "    --CHH      CHH ratio for heatmap, [0-1], default 0.1\n");

    fprintf(stderr, "    --bigwig   print bigwig files.\n");
    fprintf(stderr, "\n[align paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat align'\n");
    fprintf(stderr, "\n[calmeth paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat calmeth'\n");
    fprintf(stderr, "\n[annotion paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat annoation'\n");
    fprintf(stderr, "\n[methyPlot paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat methyPlot'\n");
    fprintf(stderr, "\n[batDMR paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat batDMR'\n");
    fprintf(stderr, "\n[visul2sample paramaters:]\n");
    fprintf(stderr, "    see the details in 'RMat visul2sample'\n\n");

    fprintf(stderr, "-h|--help   usage\n\nRMat is a naive tool, if you meet any problems or have good suggestion, please let us know. We will fix it asap!\nE-mail: qiangwei.zhou2013@gmail.com\n\n");
}

string AlignHelp = "Align help";

struct READ
{
	char Description[MAXDES];
	char Raw_Tag_Copy[MAXDES];
	char Tag_Copy[MAXTAG];
	char Quality[MAXTAG];
	char Plus[MAXTAG];
	int NCount;//Number of N's
	char N[MAXTAG];
	char NLocations[MAXTAG];
	unsigned Tag_Number;//Head =1, Tail =2
	unsigned Read_Number;
	int Real_Len;
    off64_t FLength;
};

char Read_Tag(FILE *Input_File,const char FILETYPE,READ & Read );
char Read_Tag(FILE *Input_File,FILE *Mate_File,const char FILETYPE, READ & Read,READ & Mate );
char Read_Tag_gz(gzFile Input_File,const char FILETYPE, READ & Read );
char Read_Tag_gz(gzFile Input_File,gzFile Mate_File,const char FILETYPE, READ & Read,READ & Mate);

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

void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
  std::string::size_type pos1, pos2;
  pos2 = s.find(c);
  pos1 = 0;
  while(std::string::npos != pos2)
  {
    v.push_back(s.substr(pos1, pos2-pos1));
 
    pos1 = pos2 + c.size();
    pos2 = s.find(c, pos1);
  }
  if(pos1 != s.length())
    v.push_back(s.substr(pos1));
}

// G-->A
inline void ReplaceGtoA(READ & R) //char* Read
{
	int i;
    for (i=0;i<R.Real_Len;i++)
    {
        if (R.Tag_Copy[i] == 'G' || R.Tag_Copy[i] == 'g') R.Tag_Copy[i]='A';
    }
}
//C-->T
inline void ReplaceCtoT(READ & R)
{
	int i;
    for (i=0;i<R.Real_Len;i++)
    {
        if (R.Tag_Copy[i]  == 'C' || R.Tag_Copy[i] == 'c') R.Tag_Copy[i]='T';
    }
}

//function
string get_path(string filepath);
void QCSingle(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix);
void QCPaired(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix);
void GetFileList(string PATH, FILE* outfile, string contain1);
void printoutputfiles(string outputdir, string mkpath, string output_prefix);
void FileList(string PATH, string contain1, string contain2, vector<string>& files2);
void string_replace( std::string &strBig, const std::string &strsrc, const std::string &strdst);
void fastptrim(string outputdir, string output_prefix, string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean, bool pairedend);
void runpipe(string outputdir, string output_prefix, string mkpath, string input_prefix, string input_prefix1, string input_prefix2, bool pairedend);
void get_fileformat(char* processdir, string& processname);
pthread_mutex_t meth_counter_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t flists;

static inline void trim_readno(char *s, int l)
{
    if (l > 2 && s[l-2] == '/' && isdigit(s[l-1]))
        s[l-2] = 0;
}

void convertfq(int Thread_ID, vector <string> filelist1, 
    vector <string> outfilelist1){
    int fi=0;
    gzFile gzfp_head;
    FILE *Input_head;
    READ R;
    READ R_CT, R_GA;
    char temp[MAX_PATH];
    int lentemp=0;
    string fileformat;
    int gzinfile = 0;
    for(; fi < filelist1.size(); ){
        //gzFile gzfp_head_output;
        //gzfp_head_output = gzopen(outfilelist1[fi].c_str(), "wb");
		if(Thread_ID==1){
			fprintf(stderr, "Process input file: %s\n", filelist1[fi].c_str());
		}
        gzinfile = gzfilelist[fi];
        
        if(gzinfile){
            //gzfp_head = gzopen(filelist1[fi].c_str(), "rb");
            gzfp_head = strcmp(filelist1[fi].c_str(), "-")? gzopen(filelist1[fi].c_str(), "rb") : gzdopen(fileno(stdin), "rb");
        }else{
            Input_head=File_Open(filelist1[fi].c_str(),"r");
        }
	    while ((gzinfile && Read_Tag_gz(gzfp_head,2,R) ) || (!gzinfile && Read_Tag(Input_head,2,R)))
	    {
            R.Real_Len=0;
	    	for(;R.Tag_Copy[R.Real_Len]!=0 && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
            R.Tag_Copy[R.Real_Len]='\0';
            R.Quality[R.Real_Len]='\0';
            lentemp= 0;
            for(;R.Plus[lentemp]!=0 && R.Plus[lentemp]!='\n';lentemp++);
            R.Plus[lentemp]='\0';

            R_CT=R;
            ReplaceCtoT(R_CT);
            
            fprintf(stdout, "%s YS:Z:%s YS:Z:CT\n%s\n%s\n%s\n", R_CT.Description, R.Tag_Copy, R_CT.Tag_Copy, R_CT.Plus, R_CT.Quality);
            if(ALLMAP){
                R_GA=R;
                ReplaceGtoA(R_GA);
                fprintf(stdout, "%s YS:Z:%s YS:Z:GA\n%s\n%s\n%s\n", R_GA.Description, R.Tag_Copy, R_GA.Tag_Copy, R_GA.Plus, R_GA.Quality);
            }
            //gzprintf(gzfp_head_output, "%s_YS:Z:%s_YS:Z:CT\n%s\n%s\n%s\n", R_CT.Description, R.Tag_Copy, R_CT.Tag_Copy, R_CT.Plus, R_CT.Quality);
        }
        //gzclose(gzfp_head_output);

        if(fi+1 == filelist1.size()) break;
        pthread_mutex_lock(&flists);
        fi++;
        pthread_mutex_unlock(&flists);
    }
}

void convertfq(int Thread_ID, vector <string> filelist1, vector <string> filelist2,
    vector <string> outfilelist1, vector <string> outfilelist2){
    int fi=0;
    gzFile gzfp_head, gzfp_tail;
    FILE *Input_head, *Input_tail;
    READ R, M;
    READ R_CT, M_GA;
    READ R_GA, M_CT;
    char temp[MAX_PATH];
    int lentemp=0;
    string fileformat;
    int gzinfile = 0;
    for(; fi < filelist1.size(); ){
        //gzFile gzfp_head_output, gzfp_tail_output;
        //gzfp_head_output = gzopen(outfilelist1[fi].c_str(), "wb");
        //if(PAIRED) gzfp_tail_output = gzopen(outfilelist2[fi].c_str(), "wb");
		if(Thread_ID==1){
			fprintf(stderr, "Process input file: %s, %s\n", filelist1[fi].c_str(), filelist2[fi].c_str());
		}
        
        gzinfile = gzfilelist_1[fi];
        if(gzfilelist_1[fi] != gzfilelist_2[fi]){
            fprintf(stderr, "Paired-end fastq must be in same file format, fq or gz!\n");
            exit(0);
        }

        if(gzinfile){
            gzfp_head = gzopen(filelist1[fi].c_str(), "rb");
            if(PAIRED) gzfp_tail = gzopen(filelist2[fi].c_str(), "rb");
        }else{
            Input_head=File_Open(filelist1[fi].c_str(),"r");
            if(PAIRED) Input_tail=File_Open(filelist2[fi].c_str(),"r");
        }
	    while ((gzinfile && Read_Tag_gz(gzfp_head,gzfp_tail,2,R,M) ) || (!gzinfile && Read_Tag(Input_head,Input_tail,2,R,M)))
	    {
            R.Real_Len=0;
	    	for(;R.Tag_Copy[R.Real_Len]!=0 && R.Tag_Copy[R.Real_Len]!='\n';R.Real_Len++);
            R.Tag_Copy[R.Real_Len]='\0';
            R.Quality[R.Real_Len]='\0';
            lentemp= 0;
            for(;R.Plus[lentemp]!=0 && R.Plus[lentemp]!='\n';lentemp++);
            R.Plus[lentemp]='\0';

            if(PAIRED){
                M.Real_Len=0;
                for(;M.Tag_Copy[M.Real_Len]!=0 && M.Tag_Copy[M.Real_Len]!='\n';M.Real_Len++);
                M.Tag_Copy[M.Real_Len]='\0';
                M.Quality[M.Real_Len]='\0';
                lentemp= 0;
                for(;M.Plus[lentemp]!=0 && M.Plus[lentemp]!='\n';lentemp++);
                M.Plus[lentemp]='\0';
            }

            R_CT=R;
            if(PAIRED) M_GA=M;
            ReplaceCtoT(R_CT);
            if(PAIRED) ReplaceGtoA(M_GA);
            fprintf(stdout, "%s YS:Z:%s YS:Z:CT\n%s\n%s\n%s\n", R_CT.Description, R.Tag_Copy, R_CT.Tag_Copy, R_CT.Plus, R_CT.Quality);
            if(PAIRED) fprintf(stdout, "%s YS:Z:%s YS:Z:GA\n%s\n%s\n%s\n", M_GA.Description, M.Tag_Copy, M_GA.Tag_Copy, M_GA.Plus, M_GA.Quality);
        

            if(ALLMAP){
                R_GA=R;
                if(PAIRED) M_CT=M;
                ReplaceGtoA(R_GA);
                if(PAIRED) ReplaceCtoT(M_CT);
                fprintf(stdout, "%s YS:Z:%s YS:Z:GA\n%s\n%s\n%s\n", R_GA.Description, R.Tag_Copy, R_GA.Tag_Copy, R_GA.Plus, R_GA.Quality);
                if(PAIRED) fprintf(stdout, "%s YS:Z:%s YS:Z:CT\n%s\n%s\n%s\n", M_CT.Description, M.Tag_Copy, M_CT.Tag_Copy, M_CT.Plus, M_CT.Quality);
            }

            //gzprintf(gzfp_head_output, "%s_YS:Z:%s_YS:Z:CT\n%s\n%s\n%s\n", R_CT.Description, R.Tag_Copy, R_CT.Tag_Copy, R_CT.Plus, R_CT.Quality);
            //if(PAIRED) gzprintf(gzfp_tail_output, "%s_YS:Z:%s_YS:Z:GA\n%s\n%s\n%s\n", M_GA.Description, M.Tag_Copy, M_GA.Tag_Copy, M_GA.Plus, M_GA.Quality);
        }
        //gzclose(gzfp_head_output);
        //if(PAIRED) gzclose(gzfp_tail_output);

        if(fi+1 == filelist1.size()) break;
        pthread_mutex_lock(&flists);
        fi++;
        pthread_mutex_unlock(&flists);
    }
}

std::string Split_getlast(const std::string &str, std::string sep) {
    if (str.empty()) {
        return "";
    }
    std::string ret;
    std::string temp;
    std::string::size_type begin = str.find_first_not_of(sep);
    std::string::size_type pos = 0;

    while (begin != std::string::npos) {
        pos = str.find(sep, begin);
        if (pos != std::string::npos) {
            temp = str.substr(begin, pos - begin);
            begin = pos + sep.length();
        } else {
            temp = str.substr(begin);
            begin = pos;
        }

        if (!temp.empty()) {
            ret=temp;
            temp.clear();
        }
    }
    return ret;
}

bool Split(const std::string &str, std::vector<std::string> &ret, std::string sep, int& seqlen) {
    if (str.empty()) {
        return false;
    }
    ret.clear();
    seqlen = 0;
    std::string temp;
    std::string::size_type begin = str.find_first_not_of(sep);
    std::string::size_type pos = 0;

    while (begin != std::string::npos) {
        pos = str.find(sep, begin);
        if (pos != std::string::npos) {
            temp = str.substr(begin, pos - begin);
            begin = pos + sep.length();
        } else {
            temp = str.substr(begin);
            begin = pos;
        }

        if (!temp.empty()) {
            ret.push_back(temp);
            seqlen++;
            temp.clear();
        }
    }
    return true;
}
char Char_To_CharC[256];

void Init(){
    Char_To_CharC['A']='T';Char_To_CharC['C']='G';Char_To_CharC['G']='C';Char_To_CharC['T']='A';Char_To_CharC['a']='t';Char_To_CharC['c']='g';Char_To_CharC['g']='c';Char_To_CharC['t']='a';
	Char_To_CharC['n']='n';Char_To_CharC['N']='N';
}

void processalign(){
    int maxlen = 1000;
    int buf_size = 100000;
    char myread[buf_size];
    Init();
    
    //fgets(myread, 100000, pipefile);
    //printf("%s\n", myread);
    //exit(0);
    //scanf("%s", myread); only can get the first part by '\t'
    int i = 0, j=0, k=0, readlen=0;
    //printf("== %s\n", myread);
    char head1[maxlen], head2[maxlen], head3[maxlen], newhead2[maxlen];
    char description[maxlen], rawseq[maxlen], readinf[maxlen], chrom[maxlen], newdescription[maxlen],
    matechrom[maxlen], printseq[maxlen], cigar[maxlen], seqqual[maxlen], others[maxlen];
    int pos=0, mapq=0, insertsize =0, matepos =0;
    //string myread = a;
    int flag=0;
    std::vector<std::string> seqresults;
    int seqlen=0;
    while(1){
        memset(myread, 0, buf_size);
        fgets(myread, buf_size, stdin);
        if(myread[0] == 0)
            break;
        if(myread[0]=='@'){
            sscanf(myread, "%s%s%s", head1, head2, head3);
            if(head2[0]=='S'){
                if(head2[3]=='r') continue;
                else if(head2[3]=='f') {
                    for(i=0, j=0; head2[i]!='\0' && head2[i]!='\n'; i++){
                        newhead2[j] = head2[i];
                        if(i!=3) j++;
                    }
                    newhead2[j]='\0';
                    printf("%s\t%s\t%s\n", head1, newhead2, head3);
                    continue;
                }
            }else{
                printf("%s", myread); //自带'\n'
                continue;
            }
        }else{
            sscanf(myread, "%s\t%d\t%s%d%d%s%s%d%d%*s%s", description, &flag, chrom, &pos, &mapq, cigar, matechrom, &matepos, &insertsize, seqqual);
            Split(description, seqresults, "_YS:Z:", seqlen);
            strcpy(others, Split_getlast(myread, "AS:i:").c_str());
            if(seqlen>=2){
                strcpy(newdescription, seqresults[0].c_str());
                strcpy(rawseq, seqresults[1].c_str());
                strcpy(readinf, seqresults[2].c_str());
            }
            //sscanf(description, "%s%*[_YS:Z]%s", rawseq, readinf);
            bool is_first_read = bool(flag & 0x40);
            bool is_second_read = bool(flag & 0x80);
            bool is_minus_read = bool(flag & 0x10);
            bool is_plus_read = !is_minus_read;
            bool is_mapped =  !(flag & 0x4);
            if(chrom[0]=='*' || !is_mapped) continue;
            //rm f r
            for(k=0; chrom[k]!='\0' && chrom[k]!='\n'; k++) {
                if(k>0) chrom[k-1] = chrom[k];
            }
            chrom[k-1]='\0';
            if(matechrom[0] == 'f' || matechrom[0]=='r'){
                for(k=0; matechrom[k]!='\0' && matechrom[k]!='\n'; k++) {
                    if(k>0) matechrom[k-1] = matechrom[k];
                }
                matechrom[k-1]='\0';
            }
            if(is_plus_read){
                strcpy(printseq, rawseq);
            }else{
                // reveser completment
                for(readlen=0;rawseq[readlen]!=0 && rawseq[readlen]!='\n';readlen++);
                for (i=0;i<=readlen-1;i++)
			        printseq[readlen-1-i]=Char_To_CharC[rawseq[i]];
            }
            //mappingStrand="YC:Z:GA\tYD:Z:f";
            printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tAS:i:%s", newdescription, flag, chrom, pos, mapq, cigar, 
                matechrom, matepos, insertsize, printseq, seqqual, others);
            //printf("%s\tYS:Z:---%s\tYC:Z:----%s\t%d\t%s\t%s\n", description, rawseq, readinf, flag, chrom, matechrom);
        }
    }
}

void executeCMD(const char *cmd, string outputdir, string output_prefix)
{
    char ps[1024]={0};
    FILE *ptr;
    strcpy(ps, cmd);
    fprintf(stderr, "%s\n", cmd);
    if(output_prefix != "None" && output_prefix != ""){
	    string filelogname = outputdir + output_prefix + ".run.log";
	    FILE* flog = File_Open(filelogname.c_str(), "aw");
	    fprintf(flog, "%s\n", cmd);
	    fclose(flog);
	}
    ptr=popen(ps, "w");

    if(ptr==NULL)
    {
        fprintf(stderr, "\nRun program %s error, you can run this step alone.\n", ps);
        exit(0);
    }
    pclose(ptr);
    ptr = NULL;
}

void executeCMDdir(const char *cmd, string outputdir, string output_prefix)
{
    char ps[1024]={0};
    FILE *ptr;
    strcpy(ps, cmd);
    fprintf(stderr, "%s\n", cmd);

    ptr=popen(ps, "w");

    if(ptr==NULL)
    {
        fprintf(stderr, "\nRun program %s error, you can run this step alone.\n", ps);
        exit(0);
    }
    pclose(ptr);
    ptr = NULL;
}

size_t get_executable_path( char* processdir,char* processname, size_t len)
{
	char* path_end;
	if(readlink("/proc/self/exe", processdir,len) <=0)
		return -1;
	path_end = strrchr(processdir, '/');
	if(path_end == NULL)
		return -1;
	++path_end;
	strcpy(processname, path_end);
	*path_end = '\0';
	return (size_t)(path_end - processdir);
}
string getfilename(string filelocation){
	int pos = filelocation.find_last_of('/');
	if(pos == -1) return filelocation;
	else return string(filelocation.substr(pos + 1) );
}

void printparamter1(string mkpath, string input_prefix, string input_prefix1, string input_prefix2, string outputdir, bool pairedend, string output_prefix){
	fprintf(stderr, "[ RMat ] Process Paramater file.\n");
    //paramater file
    string fparamater = mkpath + "/" + output_prefix +".Paramater.txt";
    FILE* Fparamater = File_Open(fparamater.c_str(), "w");

    string alignmode = "Single-end";
    string infiles = input_prefix;
    if(pairedend){
        alignmode = "Paired-end";
        infiles = input_prefix1 + " || " + input_prefix2;
    }
    string Sparamater = "Program\tRMat.v1\nWorkdir\t" + workdir + "\noutputdir\t" + outputdir + "\nAligner\tRMat-align\nGenome\t" + genome_index + 
    "\nOutput-prefix\t"+ output_prefix + "\nInput\t"+ infiles +  "\nAlignment-Mode\t" + alignmode;
    fprintf(Fparamater, "%s\n", Sparamater.c_str());
    fclose(Fparamater);
}

string getstring(int n)
{
	std::stringstream newstr;
	newstr<<n;
	return newstr.str();
}

string getstring(float n)
{
	std::stringstream newstr;
	newstr<<n;
	return newstr.str();
}

string configfile="";
struct SAMPLE{
	unsigned int id;
	string infile1;
	string infile2;
	string layout;
	string outprefix;
	string outdir;
	unsigned int diffid;
};
struct ARGS{
	FILE* configfp;
};
struct Threading
{
    pthread_t Thread;
    unsigned r;
    void *ret;
    ARGS Arg;
};

/*
~/methdb/Homo_sapiens/genome/memalign c2t -1 test1.fq.gz -2 test2.fq.gz -o testo -a | \
~/methdb/Homo_sapiens/genome/memalign align -g genome.fa -p 6 | \
~/methdb/Homo_sapiens/genome/memalign sam - 

~/methdb/Homo_sapiens/genome/memalign c2t -1 test1.fq.gz -2 test2.fq.gz -o testo | \
/public/home/scheng/methdb/Homo_sapiens/genome/bwame mem  -t 6 -C -p -Y /public/home/scheng/methdb/Mus_musculus_qw/genome/bwamem/mouse_reference.fa.batmeth2.fa - | samtools sort -@ 6 - -o prefix.sort.bam

** | ~/methdb/Homo_sapiens/genome/memalign sam - > test.bt2.sam 
** 改了bwamem source code, no need to modify r or f in sam, 
** but also need to process when no_directional mode.
*/
std::vector<struct SAMPLE> v_samples;
bool printbigwig = false;
string mode = "";
int main(int argc, char* argv[])
{
	string outputdir="./";
	bool pairedend=false;
	string output_prefix = "None";
	string input_prefix = "";
	string input_prefix1 = "";
	string input_prefix2 = "";
	string mkpath;
	int NTHREAD=4;
	int allthreads =24;
	bool deletelog=false;
    string genomeprefix="";

	for(int i=2;i<argc;i++)
    {
        if(!strcmp(argv[i], "-i")){
        	input_prefix= argv[++i];
        }else if(!strcmp(argv[i], "-1")){
        	input_prefix1= argv[++i];
        }else if(!strcmp(argv[i], "-2")){
        	input_prefix2= argv[++i];
            pairedend=true;
        }else if(!strcmp(argv[i], "-o"))
            output_prefix= argv[++i];
        else if(!strcmp(argv[i], "-O")){
            outputdir= argv[++i];
            if(outputdir[outputdir.length()-1] != '/')
            	outputdir+="/";
            string cmd = "mkdir -p " + outputdir;
            if(output_prefix != "None" && output_prefix != ""){
            	string rmfile = outputdir + output_prefix + ".run.log";
            	remove(rmfile.c_str());
            	deletelog=true;
            }
            if(outputdir!="./")
                executeCMDdir(cmd.c_str(), outputdir, output_prefix);
        }
        else if(!strcmp(argv[i], "-g"))
        	genome_index= argv[++i];
        else if(!strcmp(argv[i], "--gp"))
        	genomeprefix= argv[++i];
        else if(!strcmp(argv[i], "-a"))
        	ALLMAP= true;
        else if(!strcmp(argv[i], "-of"))
            outformat = argv[++i];
        else if(!strcmp(argv[i], "-p"))
            threads = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--fastp")){
            fastp = argv[++i];
            cleanreads=false;
        }else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-help")){
            usage();
            exit(0);
        }
    }

	if (argc < 2){
	    usage();
	    exit(0);
	}
    mode = argv[1];

    char workdirtmp[MAX_PATH];   
    getcwd(workdirtmp, MAX_PATH);
    workdir = workdirtmp;
    //printf("\n%s\n", workdirtmp);

    //exe location
	char processname[1024];
	char abspathtmp[MAX_PATH];
	get_executable_path(abspathtmp, processname, sizeof(abspathtmp));
	abspath= abspathtmp;
	programname = processname;

    if(mode == "align" && aligner == "bwame"){
        if(argc < 4){ 
            fprintf(stderr, "%s\n", AlignHelp.c_str());
        }
    }

    // prepare genome c2t and build index
    if(mode == "index"){
        //genome c2t
        string cmd=abspath + "genome2cg -g " + genome_index;
        if(genomeprefix == ""){
            genomeprefix = genome_index;
        }else{
            cmd = cmd + " -p " + genomeprefix;
        }
        executeCMD(cmd.c_str(), outputdir, output_prefix);
        //prepare bin and len
        cmd=abspath + "genomebinLen " + genome_index;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
        //build index
        cmd = abspath + "bwame index " + genomeprefix + ".batmeth2.fa";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
        //genome_others = genomeprefix + ".batmeth2.fa";
    }
    if(mode == "c2t"){
        // read qc and c2t / g2a
        if(input_prefix1!="" && input_prefix2!=""){
            fprintf(stderr, "Process paired-end reads!\n");
            QCPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
        if(input_prefix!=""){
            fprintf(stderr, "Process single-end reads\n");
            QCSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
    }

    if(mode == "align"){
        string cmd = abspath + "bwame mem ";
        if(ALLMAP){
            cmd += "-Z 1";
        }
        cmd = cmd + " -t " + getstring(threads) + " -C -p -Y " + genome_index + ".batmeth2.fa - ";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }

    if(mode == "sam"){
        char seq[3000];
        char seqname[100], chrom[100], cigar[100], matechrom[100];
        char readseq[300], qual[300], rawseq[300], seqconversion[3], alignSecond[100];
        char extval[1000], newextval[1000];
        int alignN = 0, alignscore = 0;
        int flag=-1, pos = 0, mapq = 0, matepos = 0, insertsize = 0, readlen = 0;
        bool is_first_read = 0;
        bool is_second_read = 0;
        bool is_minus_read = 0;
        bool is_plus_read = 0;
        bool is_mapped =  0;
        const char delim[2] = "\t";
        const char delim2[2] = ";";
        char *token = NULL;
        if(outformat == "bam"){
            samfile_t *in = 0, *out = 0;
            //bam_header_t *header;

            char out_mode[5];
            strcpy(out_mode, "wb"); //h for header
            //char *fn_in = 0, *fn_out = 0;
            char fn_in[100], fn_out[100];
            strcpy(fn_out, output_prefix.c_str());
            strcpy(fn_in, "-");

            //sam read
            //fn_in =  "-"; //(optind < argc)? argv[optind] : "-";
            // open file handlers
            if ((in = samopen(fn_in, "r", 0)) == 0) { //rb for bam file
                fprintf(stderr, "fail to open \"%s\" for reading.\n", fn_in);
            }
            
            //sam header
            if (in->header == 0) {
                fprintf(stderr, "[main_samview] fail to read the header from align file\n");
            }
            // get header from bam
            //header = bam_header_dup((const bam_header_t*) BamInFile->header);
            // modify header in->header
            /*
            header struct in bam.h
            typedef struct {
                int32_t n_targets;
                char **target_name;
                uint32_t *target_len;
                void *dict, *hash, *rg2lib;
                size_t l_text, n_text;
                char *text;
            } bam_header_t;
            */
            // write header
            if (1) { //strchr(mode, 'h')
                int i;
                // parse the header text
                fprintf(stderr, "AAAC\n");
                fprintf(stderr, "DD %d %d\n", in->header->l_text, in->header->text);

                { // then dump ->target_{name,len}
                fprintf(stderr, "AAACd2\n");
                    int j = 0;
                    for (i = 0; i < in->header->n_targets; ++i){
                        if(in->header->target_name[i][0] == 'f'){
                            strcpy(in->header->target_name[j++], in->header->target_name[i]+1);
                        }
                        //fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n", in->header->target_name[i], in->header->target_len[i]);
                    }
                    in->header->n_targets = j;
                }
                fprintf(stderr, "AAACd3 %d\n", in->header->n_targets);
            }


            //out file
            if ( (out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) {
                fprintf(stderr, "[main_samview] fail to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
                exit(0);
            }
            
            // write entire file
            bam1_t *b = bam_init1();
            int r;
            while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
                //b->core->tid
                //b->data //cigar-qname-seq-qual-aux
                samwrite(out, b); // write the alignment to `out'
            }
            if (r < -1) {
                fprintf(stderr, "[main_samview] truncated file.\n");
            }
            bam_destroy1(b);
        }else{
            // process sam and out sam
            FILE* samin = stdin;
            FILE* outsamF = File_Open(output_prefix.c_str(), "w");
            while(fscanf(samin,"%[^\n]%*c",seq)==1) //读取一行
            {
                //printf("%s\n", seq);
                if(seq[0] == '@'){
                    if(seq[1] == 'S'){
                        sscanf(seq, "@SQ\tSN:%s\tLN:%d", chrom, &pos);
                        if(chrom[0] == 'f'){
                            /*for(int j=0; j<strlen(chrom)-1; j++){
                                chrom[j] = chrom[j+1];
                            }
                            chrom[strlen(chrom)-1] = '\0';
                            */
                            strcpy(chrom,chrom+1);
                            fprintf(outsamF, "@SQ\tSN:%s\tLN:%d\n", chrom, pos);
                        }
                    }else if(seq[1] == 'P'){
                        fprintf(outsamF, "%s\n", seq);
                    }
                    continue;
                }
                sscanf(seq, "%s%d%s%d%d%s%s%d%d%s%s%[^\n]", 
                seqname, &flag, chrom, &pos, &mapq, cigar, matechrom, &matepos, &insertsize, readseq, qual, extval);
                
                if(chrom[0]=='*') {
                    fprintf(outsamF, "%s\n", seq);
                    continue;
                }
                
                is_first_read = bool(flag & 0x40);
                is_second_read = bool(flag & 0x80);
                is_minus_read = bool(flag & 0x10);
                is_plus_read = !is_minus_read;
                is_mapped =  !(flag & 0x4);
                readlen = strlen(readseq);
                
                //remove f or r in chrom
                strcpy(chrom,chrom+1);
                if(matechrom[0]!='*'){
                    strcpy(matechrom,matechrom+1); 
                }
                // split extval
                strcpy(newextval, "");
                for(token = strtok(extval, delim); token != NULL; token = strtok(NULL, delim))
                {
                    //printf("%s ", token);
                    if(token[0]=='A' && token[1] == 'S'){
                        alignscore = atoi(token+5);
                        strcat(newextval, "\t");
                        strcat(newextval, token);
                    }else if(token[0]=='Y' && token[1] == 'S'){
                        strcpy(rawseq, token+5);
                    }else if(token[0]=='Y' && token[1] == 'C'){
                        strcpy(seqconversion, token+5);
                    }else if(token[0]=='A' && token[1] == 'A'){
                        alignN = atoi(token+5);
                        strcat(newextval, "\t");
                        strcat(newextval, token);
                    }else if(token[0]=='S' && token[1] == 'A'){
                        strcpy(alignSecond, token+5);
                        strcat(newextval, "\tSA:Z:");
                        for(token = strtok(alignSecond, delim2); token != NULL; token = strtok(NULL, delim2)){
                            strcpy(token,token+1);
                            strcat(newextval, token);
                            strcat(newextval, ";");
                        }
                        //strcat(newextval, token);
                    }else{
                        strcat(newextval, "\t");
                        strcat(newextval, token);
                    }
                }
                //printf("NNNN %s\t%s\t%d\t%d\t%s\n", rawseq, seqconversion, alignscore, alignN, newextval);
                if(is_plus_read){
                    strcpy(readseq, rawseq);
                }else{
                    // reveser completment
                    for(readlen=0;rawseq[readlen]!=0 && rawseq[readlen]!='\n';readlen++);
                    for (int i=0;i<=readlen-1;i++)
                        readseq[readlen-1-i]=Char_To_CharC[rawseq[i]];
                }

                fprintf(outsamF, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n", 
                seqname, flag, chrom, pos, mapq, cigar, matechrom, matepos, insertsize, readseq, qual, newextval);
            }
            fclose(outsamF);
        }
    }
    
}

char* readline(char* s2t, int BATBUF, FILE* configFp){
	flockfile(configFp);
	char* state = fgets(s2t,BATBUF,configFp);
	funlockfile(configFp);
	return state;
}

void get_fileformat(char* processdir, string& processname)
{
	char* path_end;
	path_end = strrchr(processdir, '.');
	if(path_end == NULL)
		return;
	++path_end;
	processname= path_end;
	*path_end = '\0';
	//return (size_t)(path_end - processdir);
}

std::string get_path(string filepath)
{
	int n=filepath.find_last_of('/');
	string dirname=filepath.substr(0,n);
	//printf("\nss %s\n", dirname.c_str());
	return dirname;
}

void QCSingle(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix){
    if(aligner == "no")
        return;
    if( (input_prefix == "None") || (output_prefix == "None")){
        fprintf(stderr, "Please check the pramater.\ngenome: %s\ninput: %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix.c_str(), output_prefix.c_str());
        exit(0);
    }
    std::vector <string> infilelist;
    SplitString(input_prefix, infilelist, ",");

    string cmd="";
    //clean
    std::vector<string> cleanfilelist;
    string input_clean;
    string fileformat;
    char temp[MAX_PATH];
    int i;
    string clenfiles="";
    string outfiles=""; //for trans fq files
    string cleanname = output_prefix;
    std::vector<string> outfilelist1;
    string out1_c2t;
    for(int j=0;j<infilelist.size();j++){
        for(i=0;i<infilelist[j].length();i++)
        temp[i] = infilelist[j][i];
        temp[i]='\0';
        get_fileformat(temp, fileformat);
        if( fileformat == "gz" || fileformat == "gzip" ){
            gzfilelist[j] = 1;
            if(!cleanreads){
                input_clean = getfilename(string(temp)) + "clean.gz";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            gzfilelist[j] = 0;
            if(!cleanreads){
                input_clean = getfilename(string(temp)) + "clean.fq";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else{
            fprintf(stderr, "\n%s is a unvalid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
            exit(0);
        }
        
        string input_clean1;
        string input_clean2;
        if(!cleanreads){
            if(infilelist.size() > 0) cleanname = getfilename(input_clean);
            fprintf(stderr, "[ RMat ] raw reads: %s; clean reads: %s\n", infilelist[j].c_str(), cleanname.c_str());
            fastptrim(outputdir, output_prefix, input_prefix1, input_prefix2, infilelist[j], input_clean1, input_clean2, outputdir + cleanname, false);
            infilelist[j]=cleanname;
            cleanfilelist.push_back(cleanname);
        }
        outfilelist1.push_back(out1_c2t);
    }
    //process fq
    if(!cleanreads) convertfq(1, cleanfilelist, outfilelist1);
    else convertfq(1, infilelist, outfilelist1);

}

void QCPaired(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix){
    if(aligner == "no")
        return;
    if( (input_prefix1 == "") || (output_prefix == "") || (input_prefix2 == "")){
        fprintf(stderr, "\nError! Please check the pramater.\ngenome: %s\ninput: %s, %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix1.c_str(), input_prefix2.c_str(),output_prefix.c_str());
        exit(0);
    }
    std::vector <string> infilelist1;
    SplitString(input_prefix1, infilelist1, ",");
    std::vector <string> infilelist2;
    SplitString(input_prefix2, infilelist2, ",");
	string cmd="";
    //clean
    string input_clean1;
    string out1_c2t;
    string fileformat;
    char temp[MAX_PATH];
    int i=0;
    string clenfiles1="";
    string clenfiles2="";
    string outfiles1="";
    string outfiles2="";
    std::vector<string> cleanfilelist1;
    std::vector<string> cleanfilelist2;
    std::vector<string> outfilelist1;
    std::vector<string> outfilelist2;
    string cleanname = output_prefix;
    for(int j=0;j<infilelist1.size();j++){
        for(i=0;i<infilelist1[j].length();i++)
            temp[i] = infilelist1[j][i];
        temp[i]='\0';
        get_fileformat(temp, fileformat);
        if( fileformat == "gz" || fileformat == "gzip" ) {
            gzfilelist_1[j] = 1;
            if(!cleanreads){
                input_clean1 = getfilename(string(temp)) + "clean.gz";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            gzfilelist_1[j] = 0;
            if(!cleanreads){
                input_clean1 = getfilename(string(temp)) + "clean.fq";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else{
            fprintf(stderr, "\n%s not a valid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
            exit(0);
        }
        string input_clean;
        string input_clean2;
        
        string out2_g2a;
        for(i=0;i<infilelist2[j].length();i++)
            temp[i] = infilelist2[j][i];
        temp[i]='\0';
        get_fileformat(temp, fileformat);
        if( fileformat == "gz" || fileformat == "gzip" ) {
            gzfilelist_2[j] = 1;
            if(!cleanreads){
                input_clean2 = getfilename(string(temp)) + "clean.gz";
                out2_g2a = getfilename(string(temp)) + "clean.g2a.fq.gz";
            }else{
                out2_g2a = getfilename(string(temp)) + "g2a.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            gzfilelist_2[j] = 0;
            if(!cleanreads){
                input_clean2 = getfilename(string(temp)) + "clean.fq";
                out2_g2a = getfilename(string(temp)) + "clean.g2a.fq.gz";
            }else{
                out2_g2a = getfilename(string(temp)) + "g2a.fq.gz";
            }
        }else{
            fprintf(stderr, "\n%s not a valid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
            exit(0);
        }
        if(!cleanreads) {
            fprintf(stderr, "[ RMat ] raw reads: %s, %s; clean reads: %s, %s\n", infilelist1[j].c_str(), infilelist2[j].c_str(), input_clean1.c_str(), input_clean2.c_str());
            if(infilelist1.size() > 0) cleanname = input_clean1;
            fastptrim(outputdir, output_prefix, infilelist1[j], infilelist2[j], input_prefix, outputdir + input_clean1, outputdir + input_clean2, input_clean, true);
            infilelist1[j]=input_clean1;
            infilelist2[j]=input_clean2;
            cleanfilelist1.push_back(input_clean1);
            cleanfilelist2.push_back(input_clean2);
        }
        outfilelist1.push_back(out1_c2t);
        outfilelist2.push_back(out2_g2a);
        //alignment
        
    }
    //process fq
    if(!cleanreads) convertfq(1, cleanfilelist1, cleanfilelist2, outfilelist1, outfilelist2);
    else convertfq(1, infilelist1, infilelist2, outfilelist1, outfilelist2);

}

//{-----------------------------------  READ ROUTINES ---------------------------------------------------------
char Read_Tag_gz(gzFile Input_File,const char FILETYPE, READ & Read )
{
	pthread_mutex_lock(&Lock_gzfile);
       //flockfile(Input_File);
        if (gzgets(Input_File,Read.Description,MAXDES)!=0)// read a tag...
        {
                char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
                //gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
                if(!gzgets(Input_File,Read.Tag_Copy,MAXDES)) { printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
                if (FILETYPE == FQ)
                {
                        //gzgets(Input_File,Plus,MAXTAG);//plus
                        if(!gzgets(Input_File,Read.Plus,MAXTAG)){ printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
                        //gzgets(Input_File,Quality,MAXTAG);//phred
                        if(!gzgets(Input_File,Read.Quality,MAXTAG)){ printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
                }
                else
                {
                        Read.Quality[0]='*';Read.Quality[1]=0;
                }
                pthread_mutex_unlock(&Lock_gzfile);
        //        funlockfile(Input_File);
                return true;
        }
        else
        {
        	pthread_mutex_unlock(&Lock_gzfile);
          //      funlockfile(Input_File);
                return false;
        }
}

//{-----------------------------------  READ ROUTINES ---------------------------------------------------------
char Read_Tag(FILE *Input_File,const char FILETYPE, READ & Read )
{
        flockfile(Input_File);
        if (fgets(Read.Description,MAXDES,Input_File)!=0)// read a tag...
        {
                char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
                //gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
                if(!fgets(Read.Tag_Copy,MAXDES,Input_File)) {printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
                if (FILETYPE == FQ)
                {
                        //gzgets(Input_File,Plus,MAXTAG);//plus
                        if(!fgets(Read.Plus,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
                        //gzgets(Input_File,Quality,MAXTAG);//phred
                        if(!fgets(Read.Quality,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
                }
                else
                {
                        Read.Quality[0]='*';Read.Quality[1]=0;
                }
                funlockfile(Input_File);
                return true;
        }
        else
        {
                funlockfile(Input_File);
                return false;
        }
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_Tag
 *  Description:  Read a line from FASTA/FASTQ
 *  		  return true if successfull...
 * =====================================================================================
 */
char Read_Tag(FILE *Input_File,FILE *Mate_File,const char FILETYPE, READ & Read,READ & Mate)
{
	flockfile(Input_File);
	if(PAIRED) flockfile(Mate_File);
	Read.FLength=ftello64(Input_File);
	if (fgets(Read.Description,MAXDES,Input_File)!=0)// read a tag...
	{
		char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
		//gzgets(Input_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
		if(!fgets(Read.Tag_Copy,MAXDES,Input_File)) {printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
		if (FILETYPE == FQ)
		{
			//gzgets(Input_File,Plus,MAXTAG);//plus
			if(!fgets(Read.Plus,MAXTAG,Input_File)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
			//gzgets(Input_File,Quality,MAXTAG);//phred
			if(!fgets(Read.Quality,MAXTAG,Input_File)){ printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
		}
		else
		{
			Read.Quality[0]='*';Read.Quality[1]=0;
		}

//Mate processing...
		if (PAIRED && fgets(Mate.Description,MAXDES,Mate_File)!=0)// read a tag...
		{
			char* C=Mate.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
			//gzgets(Mate_File,Current_Tag-IGNOREHEAD,MAXDES);//tag
			if(!fgets(Mate.Tag_Copy,MAXDES,Mate_File)) { printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
			if (FILETYPE == FQ)
			{
				//gzgets(Mate_File,Plus,MAXTAG);//plus
				if(!fgets(Mate.Plus,MAXTAG,Mate_File)){ printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
				//gzgets(Mate_File,Quality,MAXTAG);//phred
				if(!fgets(Mate.Quality,MAXTAG,Mate_File)){ printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
			}
			else
			{
				Mate.Quality[0]='*';Mate.Quality[1]=0;
			}
			if(PAIRED) funlockfile(Mate_File);
		}

		funlockfile(Input_File);
		return true;
	}
	else 
	{
		if(PAIRED) funlockfile(Mate_File);
		funlockfile(Input_File);
		return false;
	}
}

char Read_Tag_gz(gzFile Input_File,gzFile Mate_File,const char FILETYPE, READ & Read,READ & Mate)
{
	//flockfile(Input_File);
	//if(PAIRED) flockfile(Mate_File);
	//Read.FLength=ftello64(Input_File);
	pthread_mutex_lock(&Lock_gzfile);
	if (gzgets(Input_File,Read.Description,MAXDES)!=0)// read a tag...
	{
		char* C=Read.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
		if(!gzgets(Input_File,Read.Tag_Copy,MAXDES)) {printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
		if (FILETYPE == FQ)
		{
			if(!gzgets(Input_File,Read.Plus,MAXTAG)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
			if(!gzgets(Input_File,Read.Quality,MAXTAG)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
		}
		else
		{
			Read.Quality[0]='*';Read.Quality[1]=0;
		}

//Mate processing...
		if (PAIRED && gzgets(Mate_File,Mate.Description,MAXDES)!=0)// read a tag...
		{
			char* C=Mate.Description;while (*C!=' ' && *C!='\t' &&*C != '\r' && *C != '\n') C++;*C=0;
			if(!gzgets(Mate_File,Mate.Tag_Copy,MAXDES)) {printf ("Read_Tag():Error reading file..\n");exit(-1);};//tag
			if (FILETYPE == FQ)
			{
				if(!gzgets(Mate_File,Mate.Plus,MAXTAG)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//plus
				if(!gzgets(Mate_File,Mate.Quality,MAXTAG)){printf ("Read_Tag():Error reading file..\n");exit(-1);};//phred
			}
			else
			{
				Mate.Quality[0]='*';Mate.Quality[1]=0;
			}
			//if(PAIRED) funlockfile(Mate_File);
		}
		pthread_mutex_unlock(&Lock_gzfile);
		//funlockfile(Input_File);
		return true;
	}
	else 
	{
		pthread_mutex_unlock(&Lock_gzfile);
		//if(PAIRED) funlockfile(Mate_File);
		//funlockfile(Input_File);
		return false;
	}
}

void GetFileList(string PATH, FILE* outfile, string contain1){
	struct dirent *ptr;
    DIR *dir;
    //string PATH = "/Users/anitafang/Downloads/Datasets/300W/300w/01_Indoor";
    dir=opendir(PATH.c_str());
    vector<string> files;
    //cout << "file list: "<< endl;
    while((ptr=readdir(dir))!=NULL)
    {
        //跳过'.'和'..'两个目录
        if(ptr->d_name[0] == '.')
            continue;
        //cout << ptr->d_name << endl;
        files.push_back(ptr->d_name);
    }
    for (int i = 0; i < files.size(); ++i)
    {

        if(outfile!=NULL)
        {
        	if(files[i].find(contain1) != string::npos)
            	fprintf(outfile, "%s\t%s\n", files[i].c_str(), PATH.c_str());    //message是程序中处理的数据
        }
        else
        {
            fprintf(stderr, "\nWarning: cat not open outfiles.!\n");
        }
    }
}

void printoutputfiles(string outputdir, string mkpath, string output_prefix){
	fprintf(stderr, "[ RMat ] Output files ...\n");
	string fouts = mkpath + "/images/outfiles.txt";
	FILE* outfiles = File_Open(fouts.c_str(), "w");
	fprintf(outfiles, "File_Name\tDirectory\n");
    GetFileList(outputdir, outfiles, output_prefix);
    string path=mkpath + "/images/";
    GetFileList(path, outfiles, "png");
    GetFileList(path, outfiles, "txt");
    GetFileList(path, outfiles, "html");
    fclose(outfiles);
}

void FileList(string PATH, string contain1, string contain2, vector<string>& files2){
	struct dirent *ptr;
    DIR *dir;
    //string PATH = "/Users/anitafang/Downloads/Datasets/300W/300w/01_Indoor";
    dir=opendir(PATH.c_str());
    vector<string> files;
    //cout << "file list: "<< endl;
    while((ptr=readdir(dir))!=NULL)
    {
        //跳过'.'和'..'两个目录
        if(ptr->d_name[0] == '.')
            continue;
        //cout << ptr->d_name << endl;
        files.push_back(ptr->d_name);
    }
    for (int i = 0; i < files.size(); ++i)
    {

        if(files[i].find(contain1) != string::npos && files[i].find(contain2) != string::npos)
        	files2.push_back(files[i]);
    }
    //return files2;
}

void string_replace( std::string &strBig, const std::string &strsrc, const std::string &strdst)
{
    std::string::size_type pos = 0;
    std::string::size_type srclen = strsrc.size();
    std::string::size_type dstlen = strdst.size();

    while( (pos=strBig.find(strsrc, pos)) != std::string::npos )
    {
        strBig.replace( pos, srclen, strdst );
        pos += dstlen;
    }
}

string string_replace_return( std::string strBig, const std::string &strsrc, const std::string &strdst)
{
    std::string::size_type pos = 0;
    std::string::size_type srclen = strsrc.size();
    std::string::size_type dstlen = strdst.size();
    string tmp=strBig;
    while( (pos=tmp.find(strsrc, pos)) != std::string::npos )
    {
        tmp.replace( pos, srclen, strdst );
        pos += dstlen;
    }
    return tmp;
}

void fastptrim(string outputdir, string output_prefix, string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean, bool pairedend){
    if(fastp != ""){
    	string cmd;
        if(pairedend)
            cmd = fastp + " -Y 0 -i " + input_prefix1 + " -I " + input_prefix2 + " -o " + input_clean1 + " -O " + input_clean2;
        else
        	cmd = fastp + " -Y 0 -i " + input_prefix + " -o " + input_clean;
        cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
}

// whole pipeline for DNA methylation analysis. Contains alignment, calute meth level, DNA methylation annatation
// on gff file or bed region, DNA methylation visulization. Differentail analysis use diffmeth function.

void runpipe(string outputdir, string output_prefix, string mkpath, string input_prefix, string input_prefix1, string input_prefix2, bool pairedend){
	fprintf(stderr, "[ RMat ] Genome: %s\n", genome_index.c_str());
	fprintf(stderr, "[ RMat ] Input file:  %s, %s %s\n", input_prefix.c_str(), input_prefix1.c_str(), input_prefix2.c_str());
	fprintf(stderr, "[ RMat ] Outfile prefix: %s\n", output_prefix.c_str());

    printparamter1(mkpath, input_prefix, input_prefix1, input_prefix2, outputdir, pairedend, output_prefix);
    fprintf(stderr, "[ RMat ] Alignment ...\n");
    if(pairedend){
        QCPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
    }
    else{
        QCSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
    }
    fprintf(stderr, "[ RMat ] Alignment summary ...\n");
    fprintf(stderr, "[ RMat ] Sorting align file ...\n");
    fprintf(stderr, "[ RMat ] ");
}
