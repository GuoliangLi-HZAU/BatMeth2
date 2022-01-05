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
#include <pthread.h>

using std::string;
using std::vector;
//alignment
int threads = 8;
float mismatches = 0.1;
string genome_index = "";
string outformat = "BAM";

string bismark2paramaters="";
string bsmapparamaters="";
string bwamethparamaters="";
// fastp
string fastp="";

//calmeth
int Qual = 30;
int redup = 1;
int region = 1000;
int sammeth = 0;
//calmeth and methyGff
int coverage = 4;
int maxcoverage = 1000;
int binCover = 1;
int chromstep = 50000;
float hs=0.4;
//methyGff
float step = 0.01;
int distance = 2000;
string gfffile = "None";
string bedfile = "None";
bool GTF = false;
//mode
string mode = "pipel";

// heatmap
float CG = 0.6;
float CHG = 0.2;
float CHH = 0.1;

//aligner
string aligner = "BatMeth2"; // bwa-meth bsmap bismark
string genome_others = "";

string abspath = "";
string workdir;
bool cleanreads=true;
string programname;

void usage(){
    fprintf(stderr, "\nBatMeth2 [mode] [paramaters]\n");
    fprintf(stderr, "mode:  index, pipel, align, calmeth, methyGff, methyPlot, batDMR, visul2sample, DMCplot\n\n");
    fprintf(stderr, "Example:\n   BatMeth2 pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -1 in1.fq.gz -2 in2.fq.gz --gff gene.gff -o outprefix\n");
    fprintf(stderr, "Or single-end:\n   BatMeth2 pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -i in.fq.gz --gff gene.gff -o outprefix\n\n");
//    ~/software_devp/batmeth2/bin/test pipel --aligner bwa-meth --go ~/practice/Genome/arabidopsis/arabidopsis_bwa-meth/TAIR10_chr_all.fa --fastp ~/software/fastp/fastp -g ~/practice/Genome/arabidopsis/arabidopsis_batmeth2_index/TAIR10_chr_all.fa -1 test1.fq -2 test2.fq --gff ~/practice/Genome/arabidopsis/TAIR10.gene.modify.gff -f 1 -o pipel.clean
    fprintf(stderr, "\n[build_index]\n");
    fprintf(stderr, "    Usage: BatMeth2 index -g genomefile. (wgbs data, must run this step first), or\n");
    fprintf(stderr, "    Usage: BatMeth2 index_rrbs -g genomefile. (rrbs data)\n");
    fprintf(stderr, "    -S    Set restriction enzyme digestion sites. default: [-S C-CGG] for MspI digestion.\n");
    fprintf(stderr, "    --gp  Genome index output prefix, default: same as genome name.\n");

    fprintf(stderr, "\n[pipel (Contains: align, calmeth, methyGff, methyPlot)]\n");
    fprintf(stderr, "[fastp location]\n");
    fprintf(stderr, "    --fastp    fastp program location.\n");
    fprintf(stderr, "               ** If --fastp is not defined, the input file should be clean data.\n");
    //fprintf(stderr, "\n[select aligner]\n");
    //fprintf(stderr, "    --aligner         BatMeth2(default), bwa-meth(v1), bsmap, bismark2, \n");
    //fprintf(stderr, "                      or no (exit output_prefix.sam file, no need align again)\n");
    //fprintf(stderr, "    --bismark2para    bismark2 paramaters, contained by \"paramaters\"\n");
    //fprintf(stderr, "    --bsmsppara       bsmap paramaters\n");
    //fprintf(stderr, "    --bwamethpara     bwameth paramaters\n");
    //fprintf(stderr, "[other aligners paramaters]\n");
    //fprintf(stderr, "    --go              Name of the genome, contaion index build by aligner. (bwa-meth/bismark2)\n");
    fprintf(stderr, "\n[main paramaters]\n");
    //fprintf(stderr, "    --config [config file].   When we run pipel function in batches datasets, \n");
    //fprintf(stderr, "                              please fill in the specified configuration file. \n");
    //fprintf(stderr, "                              And there is a sample file (multirun.onf) in the BatMeth2 directory.\n");
    //fprintf(stderr, "    --mp [4]                  When batch processing data, we set the number of samples to run at a time (-mp, default is 4), and each sample needs six threads (- P parameter) by default.\n");
    //fprintf(stderr, "    --ap [24]                 When batch processing data, the number of samples to run at a time (default 4) * threads of each sample (default 6)\n");
    fprintf(stderr, "    -o [outprefix]            Name of output file prefix\n");
    fprintf(stderr, "    -O [out folder]           Output of result file to specified folder, default output to current folder (./)\n");
    fprintf(stderr, "    -of [SAM/BAM]             Output format, default BAM.\n");
    fprintf(stderr, "\n[alignment paramaters]\n");
    fprintf(stderr, "    -i    Name of input file, support .fq/.fastq and .gz/.gzip format. if paired-end. please use -1, -2\n");
    fprintf(stderr, "    -1    Name of input file left end, if single-end. please use -i\n");
    fprintf(stderr, "    -2    Name of input file right end\n");
    fprintf(stderr, "          -i/-1/-2 can be comma-separated lists (no whitespace), only supported in BatMeth2 aligner.\n");
    fprintf(stderr, "    -g    Name of the genome mapped against, make sure build index first.\n");
    fprintf(stderr, "    -p    Launch <integer> threads\n");
    fprintf(stderr, "\n[calmeth paramaters]\n");
    fprintf(stderr, "    --Qual      calculate the methratio while read QulityScore >= Q. default:30\n");
    fprintf(stderr, "    --redup     REMOVE_DUP, 0 or 1, default 1\n");
    fprintf(stderr, "    --region    Bins for DMR calculate , default 1000bp .\n");
    fprintf(stderr, "    -f          for sam format outfile contain methState. [0 or 1], default: 0 (dont output this file).\n");
	fprintf(stderr, "    -n          maximum mismatches allowed due to seq. default 0.1 percentage of the read length. [0-0.3]\n");
    fprintf(stderr, "\n[calmeth and methyGff paramaters]\n");
    fprintf(stderr, "    --coverage    >= <INT> coverage. default:4\n");
    fprintf(stderr, "    --binCover    >= <INT> nCs per region. default:1\n");
    fprintf(stderr, "    --chromstep   Chromsome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp)\n");
    fprintf(stderr, "\n[methyGff paramaters]\n");
    fprintf(stderr, "    --gtf/--gff/--bed    Gtf or gff file / bed file\n");
    fprintf(stderr, "    --distance           DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000\n");
    fprintf(stderr, "    --step               Gene body and their flanking sequences using an overlapping sliding window of 0.02 of the sequence length at a step of 0.01 of the sequence length. So default step: 0.01 (1%) \n");
    fprintf(stderr, "    -C                   <= <INT> coverage. default:1000\n");
    fprintf(stderr, "\n[MethyPlot paramaters]\n");
    fprintf(stderr, "    --CG       CG ratio for heatmap, [0-1], default 0.6\n");
    fprintf(stderr, "    --CHG      CHG ratio for heatmap, [0-1], default 0.2\n");
    fprintf(stderr, "    --CHH      CHH ratio for heatmap, [0-1], default 0.1\n");

    fprintf(stderr, "    --bigwig   print bigwig files.\n");
    fprintf(stderr, "\n[align paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 align'\n");
    fprintf(stderr, "\n[calmeth paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 calmeth'\n");
    fprintf(stderr, "\n[methyGff paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 methyGff'\n");
    fprintf(stderr, "\n[methyPlot paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 methyPlot'\n");
    fprintf(stderr, "\n[batDMR paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 batDMR'\n");

    fprintf(stderr, "-h|--help   usage\n\nBatMeth2 is a naive tool, if you meet any problems or have good suggestion, please let us know. We will fix it asap!\nE-mail: qiangwei.zhou2013@gmail.com\n\n");
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

void onlyexecuteCMD(const char *cmd, string errorinfor)
{
    char ps[1024]={0};
    FILE *ptr;
    strcpy(ps, cmd);
    //fprintf(stderr, "[MM] %s\n", cmd);
    ptr=popen(ps, "w");

    if(ptr==NULL)
    {
        fprintf(stderr, "\n%s\n", errorinfor.c_str());
        exit(0);
    }
    pclose(ptr);
    ptr = NULL;
}

void executeCMD(const char *cmd, string outputdir, string output_prefix)
{
    char ps[1024]={0};
    FILE *ptr;
    strcpy(ps, cmd);
    fprintf(stderr, "[MM] %s\n", cmd);
    if(output_prefix != "None" && output_prefix != ""){
	    string filelogname = outputdir + output_prefix + ".cmd.log";
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
	fprintf(stderr, "[ BatMeth2 ] Process Paramater file.\n");
    //paramater file
    string fparamater = mkpath + "/" + output_prefix +".Paramater.txt";
    FILE* Fparamater = File_Open(fparamater.c_str(), "w");

    string alignmode = "Single-end";
    string infiles = input_prefix;
    if(pairedend){
        alignmode = "Paired-end";
        infiles = input_prefix1 + " || " + input_prefix2;
    }
    string Sparamater = "Program\tBatMeth2.v1\nWorkdir\t" + workdir + "\noutputdir\t" + outputdir + "\nAligner\tBatMeth2-align\nGenome\t" + genome_index + 
    "\nAnnotation\t" + gfffile + "/" + bedfile + "\nOutput-prefix\t"+ output_prefix + "\nInput\t"+ infiles +  "\nAlignment-Mode\t" + alignmode;
    fprintf(Fparamater, "%s\n", Sparamater.c_str());
    fclose(Fparamater);
}

void printparamter2(string mkpath, string output_prefix){
	fprintf(stderr, "[ BatMeth2 ] Process Paramater file2.\n");
    //paramater2
    string fparamater2 = mkpath + "/" + output_prefix + ".Paramater2.txt";
    FILE* Fparamater = File_Open(fparamater2.c_str(), "w");
    fprintf(Fparamater, "Threads\t%d\nCalmeth\t----\nQuality_Score\t%d\nredup\t%d\nmeth region length\t%d\nPrint methstate samfile\t%d\ncalmeth and methyGff\t----\nCoverage\t%d\nmaxCoverage\t%d\nbinCoverage\t%d\nchromStep\t%d\nmethyGff\t----\nGene bins step\t%.3f\nDistance of upstream and downstream\t%d  ", threads, Qual, redup, region, sammeth, coverage, maxcoverage, binCover, chromstep, step, distance);
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

//function
void calmeth(string inputf, string outputdir, string output_prefix);
void build_index(string para, string command);
string get_path(string filepath);
void alignmentSingle(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix);
void alignmentPaired(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix);
void annotation(string outputdir, string output_prefix);
void GetFileList(string PATH, FILE* outfile, string contain1);
void printoutputfiles(string outputdir, string mkpath, string output_prefix);
void FileList(string PATH, string contain1, string contain2, vector<string>& files2);
void string_replace( std::string &strBig, const std::string &strsrc, const std::string &strdst);
void mvpng(string outputdir, string mkpath, string output_prefix);
void doc2html(string outputdir, string mkpath, string output_prefix);
void alignmentstate(string outputdir, string output_prefix, string mkpath);
void methyPlot(string outputdir, string output_prefix);
void fastptrim(string outputdir, string output_prefix, string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean, bool pairedend);
void runpipe(string outputdir, string output_prefix, string mkpath, string input_prefix, string input_prefix1, string input_prefix2, bool pairedend);
void detect_mode(string mode, int Nparas, char* paramaters[], string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix, string mkpath, bool pairedend, \
  string redss, string genome_index, string genomeprefix);
void get_fileformat(char* processdir, string& processname);
int MAX_PATH = 1000;
void *nprunpipel(void *arg);
pthread_mutex_t meth_counter_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

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

string bt2modeHelp = "BatMeth2 \n" \
                     "index          prepare genome index, before analysis!\n" \
                     "pipel          pipeline of DNA methylation data analysis \n" \
                     "align          data alignment\n" \
                     "calmeth        calculate DNA methylation level based on alignment file\n" \
                     "methyGff       calculate DNA methylation distribution across gene/TE/peak/other regions based on methratio file\n" \
                     "batDMR         calculate differential DNA methylation cytosine or region (DMC/DMR)\n" \
                     "methyPlot      visulization of DNA methylation data\n" ;

void indexUsage(){

    fprintf(stderr, "BatMeth2 index\n");
    fprintf(stderr, "    Usage: BatMeth2 index -g genomefile. (wgbs data, must run this step first), or\n");
    fprintf(stderr, "    Usage: BatMeth2 index_rrbs -g genomefile. (rrbs data)\n");
    fprintf(stderr, "    -S    Set restriction enzyme digestion sites. default: [-S C-CGG] for MspI digestion.\n");
    fprintf(stderr, "    --gp  Genome index output prefix, default: same as genome name.\n");
}

std::vector<struct SAMPLE> v_samples;
bool printbigwig = false;
int main(int argc, char* argv[])
{
	string outputdir="./";
	bool pairedend=false;
	string output_prefix = "None";
	string input_prefix = "None";
	string input_prefix1 = "";
	string input_prefix2 = "";
    string genomeprefix = "";
	string mkpath;
	int NTHREAD=4;
	int allthreads =24;
	bool deletelog=false;
    string redss = "C-CGG";

	for(int i=1;i<argc;i++)
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
        else if(!strcmp(argv[i], "--config"))
            configfile= argv[++i];
        else if(!strcmp(argv[i], "-O")){
            outputdir= argv[++i];
            if(outputdir[outputdir.length()-1] != '/')
            	outputdir+="/";
            string cmd = "mkdir -p " + outputdir;
            if(output_prefix != "None" && output_prefix != ""){
            	string rmfile = outputdir + output_prefix + ".run.log";
            	remove(rmfile.c_str());
            	deletelog=true;
                rmfile = outputdir + output_prefix + ".cmd.log";
            	remove(rmfile.c_str());
            	deletelog=true;
            }
            if(outputdir!="./" && outputdir!=".")
                executeCMDdir(cmd.c_str(), outputdir, output_prefix);
        }
        else if(!strcmp(argv[i], "-g"))
        	genome_index= argv[++i];
        else if(!strcmp(argv[i], "--gp"))
        	genomeprefix= argv[++i];
        else if(!strcmp(argv[i], "-of"))
                outformat = argv[++i];
        else if(!strcmp(argv[i], "-p"))
            threads = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--mp"))
            NTHREAD = atoi(argv[++i]);
        else if(!strcmp(argv[i], "-n")){
            mismatches = atof(argv[++i]);
        }
        else if(!strcmp(argv[i], "-f"))
            sammeth = atoi(argv[++i]);
        else if(!strcmp(argv[i], "-C"))
            maxcoverage = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--gtf")){
            gfffile = argv[++i];
            GTF = true;
        }else if(!strcmp(argv[i], "--fastp")){
            fastp = argv[++i];
            cleanreads=false;
        }else if(!strcmp(argv[i], "--bismark2para")){
            for(int i=1;i<strlen(argv[++i])-1;i++){
            	bismark2paramaters[i-1] = argv[++i][i];
            }
        }else if(!strcmp(argv[i], "--bsmappara")){
            for(int i=1;i<strlen(argv[++i])-1;i++){
            	bsmapparamaters[i-1] = argv[++i][i];
            }
        }else if(!strcmp(argv[i], "--bwamethpara")){
            for(int i=1;i<strlen(argv[++i])-1;i++){
            	bwamethparamaters[i-1] = argv[++i][i];
            }
        }
        else if(!strcmp(argv[i], "--gff"))
            gfffile = argv[++i];
        else if(!strcmp(argv[i], "--bed"))
            bedfile = argv[++i];
        else if(!strcmp(argv[i], "--coverage"))
            coverage = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--binCover"))
            binCover = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--chromstep"))
            chromstep = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--step"))
            step = atof(argv[++i]);
        else if(!strcmp(argv[i], "--distance"))
            distance = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--Qual"))
            Qual = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--CG"))
            CG = atof(argv[++i]);
        else if(!strcmp(argv[i], "--CHG"))
            CHG = atof(argv[++i]);
        else if(!strcmp(argv[i], "--CHH"))
            CHH = atof(argv[++i]);
        else if(!strcmp(argv[i], "--region"))
            region = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--redup"))
            redup = atoi(argv[++i]);
        else if(!strcmp(argv[i], "--aligner"))
            aligner = argv[++i];
        else if(!strcmp(argv[i], "--go"))
            genome_others = argv[++i];
        else if(!strcmp(argv[i], "--bigwig"))
            printbigwig = true;
        else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-help")){
            usage();
            exit(0);
        }else if(!strcmp(argv[i], "-S")){
			redss = argv[++i];
		}
    }

	if (argc < 2){
	    usage();
	    exit(0);
	}

    mode = argv[1];
    if(mode != "index" && mode != "index_rrbs" && mode != "pipel" && mode != "align" && mode != "calmeth" && mode != "methyGff" && mode != "methyPlot" && mode != "batDMR" && 
    mode != "visul2sample" && mode != "DMCplot"){
    	fprintf(stderr, "\nNot a valid mode\n");
    	usage();
    	exit(0);
    }
    if(fastp!=""){
        string cmd = "which " + fastp;
        onlyexecuteCMD(cmd.c_str(), "Undetected fastp command, please check the location of fastp!!");
        //exit(0); 
    }
    // workdir location
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

	if(!deletelog && output_prefix != "None" && output_prefix != ""){
		string rmfile = outputdir + output_prefix + ".run.log";
		remove(rmfile.c_str());
		deletelog=true;
        rmfile = outputdir + output_prefix + ".cmd.log";
		remove(rmfile.c_str());
	}

    if(mode == "align" && aligner == "BatMeth2"){
        if(argc < 4){ 
            string cmd = abspath + "batmeth2-align";
            executeCMD(cmd.c_str(), outputdir, output_prefix);
            exit(0);
        }
    }
    //for ann and bin file
    if(mode == "pipel" ){
        if(genome_index == "" && genome_others != ""){
        	if(aligner == "BatMeth2"){
        		fprintf(stderr, "\nPlease defined genome_index location, when use batmeth2 aligner\n");
        		exit(0);
        	}
        	string cmd = abspath + "preGenome " + genome_others;
        	executeCMD(cmd.c_str(), outputdir, output_prefix);
        	genome_index = genome_others;
        }
        if(gfffile == "None" && bedfile == "None"){
            fprintf(stderr, "\nPlease defined gff/gtf/bed file for methyGff analysis!\n");
            exit(0);
        }
	}
    mkpath = outputdir;

    if(mode == "mode"){
        fprintf(stderr, "%s\n", bt2modeHelp.c_str());
    }
    
	if (argc < 4){
	    if (mode == "pipel")
	        usage();
	    else
	        detect_mode(mode, argc, argv, outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, mkpath, pairedend, 
                redss, genome_index, genomeprefix);
	    exit(0);
	}

    if( mode != "index" && mode != "index_rrbs" ){
        fprintf(stderr, "[ Program directory ] %s\n[ Program name ] %s\n",abspathtmp, processname);
        fprintf(stderr, "[ Workdir ] %s\n", workdirtmp);
        fprintf(stderr, "[ outputdir ] %s\n", outputdir.c_str());
    }

	if(configfile!="" && mode != "pipel"){
		fprintf(stderr, "\nConfigure files only valid in pipel mode!\n");
		exit(0);
	}

	if(configfile != ""){

		fprintf(stderr, "[ BatMeth2 ] Genome: %s\n", genome_index.c_str());
		fprintf(stderr, "[ BatMeth2 ] Annotation, gtf: %s; bed: %s;\n", gfffile.c_str(), bedfile.c_str());

		FILE* configFp = File_Open(configfile.c_str(), "r");
		ARGS arg; arg.configfp = configFp;

		Threading* Thread_Info=(Threading*) malloc(sizeof(Threading)*NTHREAD);
	    pthread_attr_t Attrib;
	    pthread_attr_init(&Attrib);
	    pthread_attr_setdetachstate(&Attrib, PTHREAD_CREATE_JOINABLE);

	    for (int i=0;i<NTHREAD;i++)
	    {
	        Thread_Info[i].Arg=arg;
	        Thread_Info[i].r=pthread_create(&Thread_Info[i].Thread,&Attrib, nprunpipel ,(void*) &Thread_Info[i].Arg);
	        if(Thread_Info[i].r) {printf("Launch_Threads():Cannot create thread..\n");exit(-1);}
	    }
	    pthread_attr_destroy(&Attrib);
	    for (int i=0;i<NTHREAD;i++)
	    {
	        pthread_join(Thread_Info[i].Thread,NULL);
	    }
	    free(Thread_Info);		
		fclose(configFp);

	}else
		detect_mode(mode, argc, argv, outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, mkpath, pairedend, \
            redss, genome_index, genomeprefix);

}
char* readline(char* s2t, int BATBUF, FILE* configFp){
	flockfile(configFp);
	char* state = fgets(s2t,BATBUF,configFp);
	funlockfile(configFp);
	return state;
}

void *nprunpipel(void *arg){
	ARGS *args = (struct ARGS*)arg;
	FILE* configFp = args->configfp;
	struct SAMPLE sample; bool pairedend=true;
	int BATBUF=2000;char s2t[BATBUF];
	unsigned int sampleid;
	char layout[100];
	char infile1[2000];
	char infile2[2000];
	char outp[1000];
	char outd[1000];
	unsigned int diffsample;
	while ( readline(s2t, BATBUF, configFp)!=0 ){
			if(s2t[0] == '#') continue;
			sscanf(s2t, "%u%s%s%s%s%s%u", &sampleid, layout, infile1, infile2, outp, outd);
			if(outd[strlen(outd)-1] != '/'){
				outd[strlen(outd)]='/'; 
				outd[strlen(outd)+1]='\0';
			}
			string cmd = "mkdir -p " + (string)outd;
            executeCMD(cmd.c_str(), outd, outp);

			sample.id = sampleid;
			sample.layout = layout;
			sample.infile1 = infile1;
			sample.infile2 = infile2;
			sample.outprefix = outp;
			sample.outdir = outd;
			//sample.diffid = diffsample;
			v_samples.push_back(sample);
			if(sample.layout == "PE") pairedend = true;
			else pairedend=false;
			fprintf(stderr, "[ BatMeth2 ] Input file:  %s %s\n", infile1, infile2);
			fprintf(stderr, "[ BatMeth2 ] Outfile prefix: %s\n", outp);
                        string mkpath = (string)outd;
			runpipe(outd, outp, mkpath, infile1, infile1, infile2, pairedend);
	}
}

void *nprunpipel2(void *arg){
    SAMPLE *args = (struct SAMPLE*)arg;
    string input_prefix = args->infile1;
    string input_prefix1 = args->infile1;
    string input_prefix2 = args->infile2;
    string outputdir = args->outdir;
    string output_prefix = args->outprefix;
    bool pairedend = true;
    if(args->layout == "PE") pairedend = true;
    else pairedend = false;

    fprintf(stderr, "[ BatMeth2 ] Input file:  %s, %s %s\n", input_prefix.c_str(), input_prefix1.c_str(), input_prefix2.c_str());
    fprintf(stderr, "[ BatMeth2 ] Outfile prefix: %s\n", output_prefix.c_str());
    //string mkpath = outputdir + "/batmeth2_report_" + output_prefix;
    string mkpath = outputdir;

    printparamter1(mkpath, input_prefix, input_prefix1, input_prefix2, outputdir, pairedend, output_prefix);
    printparamter2(mkpath, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Alignment ...\n");
    string align_result = outputdir + output_prefix + ".sam";
    if(pairedend){
        alignmentPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
    }
    else{
        alignmentSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
    }
    fprintf(stderr, "[ BatMeth2 ] Alignment summary ...\n");
    fprintf(stderr, "[ BatMeth2 ] Sorting align file ...\n");
    fprintf(stderr, "[ BatMeth2 ] ");
    alignmentstate(outputdir, output_prefix, mkpath);
    fprintf(stderr, "[ BatMeth2 ] Methylation Status ...\n");
    align_result = outputdir + output_prefix + ".sort.bam";
    calmeth(align_result, outputdir, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Annotation ...\n");
    annotation(outputdir, output_prefix);
    string methratioLogfile = outputdir + output_prefix + ".log.txt";
    string newlogfile = mkpath + output_prefix + ".methbasci.txt";
    string cmd = "cp ";
    cmd += methratioLogfile; cmd+=" ";
    cmd += newlogfile;
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Visulization ...\n");
    methyPlot(outputdir, output_prefix );
    //mvpng(outputdir, mkpath, output_prefix);
    //printoutputfiles(outputdir, mkpath, output_prefix);
    //doc2html(outputdir, mkpath, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Done!\nBatMeth2 is a naive tools, if you meet any problems, please let us know. We will fix it asap!\nE-mail: qiangwei.zhou2013@gmail.com\n");
}

void calmeth(string inputf, string outputdir, string output_prefix){
    if(genome_index == ""){
        fprintf(stderr, "Please use -g pramater to defined the genome index location.\n");
        exit(0);
    }
    if(inputf == "None"){
        fprintf(stderr, "Must have input sam file (-i paramater).\n");
        exit(0);
    }
    string rd = " ";
    if(redup == 1){
        rd = " -r";
    }
    string sammethfile = " ";
    if(sammeth == 1)
        sammethfile =  " -f " + outputdir + output_prefix + ".md.sam";
    
    string cmd = abspath + "calmeth" + " -g " + genome_index + " -n " + getstring(mismatches) + " -b " + inputf + " -m " + outputdir + output_prefix  + rd + sammethfile;
    if(chromstep != 50000)
        cmd = cmd + " -s " + getstring(chromstep);
    if(coverage != 4)
        cmd = cmd + " -c " + getstring(coverage);
    if(binCover != 1)
        cmd = cmd + " -nC " + getstring(binCover);
    if(Qual != 20)
        cmd = cmd + " -Q " + getstring(Qual);
    if(region != 1000)
        cmd = cmd + " -R " + getstring(region);
    cmd = cmd + " -as 1 ";
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    return;
}

void build_index(string para, string mode, string outputdir, string output_prefix, string redss, \
    string genome_index, string genomeprefix){
    if(genome_index==""){
        indexUsage();
        exit(0);
    }
    string cmd="";
    //rrbs
    if(mode == "index_rrbs"){
        cmd = abspath + "build_index_rrbs -S " + redss + " " + genome_index;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
        genome_index = genome_index + ".rrbs.fa";
        //genome.rrbs.fa
    }
    //genome c2t
    cmd=abspath + "genome2cg -g " + genome_index;
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
    return;

    cmd = abspath + "memalign " + mode + " " + para; //sys.argv[2]
    string temp = "None";
    executeCMD(cmd.c_str(), temp, temp);
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

void alignmentSingle(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix){
    if(aligner == "no")
        return;
    if((genome_index == "" && genome_others == "") || (input_prefix == "None") || (output_prefix == "None")){
        fprintf(stderr, "Please check the pramater.\ngenome: %s\ninput: %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix.c_str(), output_prefix.c_str());
        if(aligner == "BatMeth2"){
        	string cmd = abspath + programname;
            executeCMD(cmd.c_str(), outputdir, output_prefix);
        }
        exit(0);
    }
    std::vector <string> infilelist;
    SplitString(input_prefix, infilelist, ",");

    string cmd="";
	if(!cleanreads){
		std::vector<string> cleanfilelist;
    	string input_clean;
    	string fileformat;
    	char temp[MAX_PATH];
    	int i;
    	string clenfiles="";
    	string cleanname = output_prefix;
    	for(int j=0;j<infilelist.size();j++){
	    	for(i=0;i<infilelist[j].length();i++)
			temp[i] = infilelist[j][i];
			temp[i]='\0';
	    	get_fileformat(temp, fileformat);
	    	if( fileformat == "gz" || fileformat == "gzip" ) 
	    		input_clean = string(temp) + "clean.gz";
	    	else if(fileformat == "fq" || fileformat == "fastq")
	    		input_clean = string(temp) + "clean.fq";
	    	else{
	    		fprintf(stderr, "\n%s is a unvalid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
	    		exit(0);
	    	}
	    	if(infilelist.size() > 1)
				cleanname = getfilename(input_clean);
	    	string input_clean1;
	    	string input_clean2;
	    	fprintf(stderr, "[ BatMeth2 ] raw reads: %s; clean reads: %s\n", infilelist[j].c_str(), cleanname.c_str());
		    fastptrim(outputdir, output_prefix, input_prefix1, input_prefix2, infilelist[j], input_clean1, input_clean2, outputdir + cleanname, false);
			infilelist[j]=cleanname;
			cleanfilelist.push_back(cleanname);
			//alignment
			fprintf(stderr, "[ BatMeth2 ] Alignment %s ...\n", input_clean.c_str());
		    if(aligner == "BatMeth2")
		    	if(clenfiles == "")
		    		clenfiles = outputdir + cleanname;
		    	else clenfiles = clenfiles + "," + outputdir + cleanname;
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + outputdir + cleanname + " -t " + getstring(threads) + " --prefix " + outputdir + cleanname + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + outputdir + cleanname + " -d " + genome_index + " -o " + outputdir + cleanname + ".sam -p " + getstring(threads);// + " -r 0";
		    else if(aligner == "bismark2") //need test get_path
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters + " -U " + outputdir + cleanname + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + cleanname + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    executeCMD(cmd.c_str(), outputdir, output_prefix);
	    }
	    if(aligner=="BatMeth2"){
                if(outformat == "BAM")
	    	    cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles + " | samtools sort -m 2G -@ " + getstring(threads) + " -O BAM -o " + outputdir + output_prefix + ".sort.bam -";
                else
                    cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles + " -o " + outputdir + output_prefix + ".sam";
	    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
	    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    }else{
	    	if(infilelist.size() > 1) {
		    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
		    	for(int j=0; j< cleanfilelist.size(); j++){
		    		cmd = cmd + " " + outputdir + cleanfilelist[j] + ".sam";
		    	}
		    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    	}
	    }
	}else{ // already clean data
		string rawname = output_prefix;
		for(int j=0;j<infilelist.size();j++){
			if(infilelist.size() > 1)
				rawname = getfilename(infilelist[j]);
			fprintf(stderr, "[ BatMeth2 ] Alignment %s ...\n", infilelist[j].c_str());
		    if(aligner == "BatMeth2")
		    	cmd="";
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + infilelist[j] + " -t " + getstring(threads) + " --prefix " + outputdir + rawname + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + infilelist[j] + " -d " + genome_index + " -o " + outputdir + rawname + ".sam -p " + getstring(threads);// + " -r 0";
		    else if(aligner == "bismark2") //need test get_path
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters + " -U " + infilelist[j] + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + rawname + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    executeCMD(cmd.c_str(), outputdir, output_prefix);
		}
	    if(aligner=="BatMeth2"){
                if(outformat == "BAM")
	    	    cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix + " | samtools sort -m 2G -@ " + getstring(threads) + " -O BAM -o " + outputdir + output_prefix + ".sort.bam -";
                else
                    cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix + " -o " + outputdir + output_prefix + ".sam";
	    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
	    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    }else{
	    	if(infilelist.size() > 1) {
		    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
		    	for(int j=0; j< infilelist.size(); j++){
		    		cmd = cmd + " " + outputdir + getfilename(infilelist[j]) + ".sam";
		    	}
		    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    	executeCMD(cmd.c_str(), outputdir, output_prefix);
		    	string rmfile;
		    	//for(int j=0; j< infilelist.size(); j++){
		    	//	rmfile= getfilename(infilelist[j]) + ".sam";
		    	//	remove(rmfile)
		    	//}
	    	}
	    }
	}


}

void alignmentPaired(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix){
    if(aligner == "no")
        return;
    if((genome_index == "" && genome_others == "") || (input_prefix1 == "") || (output_prefix == "None") || (input_prefix2 == "")){
        fprintf(stderr, "\nError! Please check the pramater.\ngenome: %s\ninput: %s, %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix1.c_str(), input_prefix2.c_str(),output_prefix.c_str());
        if(aligner == "BatMeth2"){
        	string cmd = abspath + programname;
            executeCMD(cmd.c_str(), outputdir, output_prefix);
        }
        exit(0);
    }
    std::vector <string> infilelist1;
    SplitString(input_prefix1, infilelist1, ",");
    std::vector <string> infilelist2;
    SplitString(input_prefix2, infilelist2, ",");
	string cmd="";
    if(!cleanreads){
    	string input_clean1;
    	string fileformat;
    	char temp[MAX_PATH];
    	int i=0;
    	string clenfiles1="";
    	string clenfiles2="";
    	std::vector<string> cleanfilelist1;
    	std::vector<string> cleanfilelist2;
    	string cleanname = output_prefix;
    	for(int j=0;j<infilelist1.size();j++){
	    	for(i=0;i<infilelist1[j].length();i++)
				temp[i] = infilelist1[j][i];
			temp[i]='\0';
	    	get_fileformat(temp, fileformat);
	    	if( fileformat == "gz" || fileformat == "gzip" ) 
	    		input_clean1 = getfilename(string(temp)) + "clean.gz";
	    	else if(fileformat == "fq" || fileformat == "fastq")
	    		input_clean1 = getfilename(string(temp)) + "clean.fq";
	    	else{
	    		fprintf(stderr, "\n%s not a valid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
	    		exit(0);
	    	}
	    	string input_clean;
	    	string input_clean2;
	    	for(i=0;i<infilelist2[j].length();i++)
				temp[i] = infilelist2[j][i];
			temp[i]='\0';
	    	get_fileformat(temp, fileformat);
	    	if( fileformat == "gz" || fileformat == "gzip" ) 
	    		input_clean2 = getfilename(string(temp)) + "clean.gz";
	    	else if(fileformat == "fq" || fileformat == "fastq")
	    		input_clean2 = getfilename(string(temp)) + "clean.fq";
	    	else{
	    		fprintf(stderr, "\n%s not a valid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
	    		exit(0);
	    	}
	    	fprintf(stderr, "[ BatMeth2 ] raw reads: %s, %s; clean reads: %s, %s\n", infilelist1[j].c_str(), infilelist2[j].c_str(), input_clean1.c_str(), input_clean2.c_str());
		    if(infilelist1.size() >1)
				cleanname = input_clean1;

		    fastptrim(outputdir, output_prefix, infilelist1[j], infilelist2[j], input_prefix, outputdir + input_clean1, outputdir + input_clean2, input_clean, true);
			infilelist1[j]=input_clean1;
			infilelist2[j]=input_clean2;
			cleanfilelist1.push_back(input_clean1);
			cleanfilelist2.push_back(input_clean2);
			//alignment
			fprintf(stderr, "[ BatMeth2 ] Alignment %s, %s...\n", input_clean1.c_str(), input_clean2.c_str());
		    if(aligner == "BatMeth2"){
		    	if(clenfiles1 == "")
		    		clenfiles1 = outputdir + input_clean1;
		    	else clenfiles1 = clenfiles1 + "," + outputdir + input_clean1;
		    	if(clenfiles2 == "")
		    		clenfiles2 = outputdir + input_clean2;
		    	else clenfiles2 = clenfiles2 + "," + outputdir + input_clean2;
			}
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + outputdir + input_clean1 + " " + outputdir + input_clean2 + " -t " + getstring(threads) + " --prefix " + outputdir + cleanname + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + outputdir + input_clean1 + " -b " + outputdir + input_clean2 + " -d " + genome_index + " -o " + outputdir + cleanname + ".sam -p " + getstring(threads); //+ " -r  0";
		    else if (aligner == "bismark2")
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters  + " -1 " + outputdir + input_clean1 + " -2 " + outputdir + input_clean2 + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + cleanname + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    executeCMD(cmd.c_str(), outputdir, output_prefix);
		}
		if(aligner=="BatMeth2"){
            if(outformat=="BAM")
                cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles1 + " -i " + clenfiles2 + " | samtools sort -m 2G -@ " + getstring(threads) + " -O BAM -o " + outputdir + output_prefix + ".sort.bam -";
            else
                cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles1 + " -i " + clenfiles2 + " -o " + outputdir + output_prefix + ".sam";
	    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
	    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    }else{
	    	if(infilelist1.size() > 1){
		    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
		    	for(int j=0; j< cleanfilelist1.size(); j++){
		    		cmd = cmd + " " + outputdir + cleanfilelist1[j] + ".sam";
		    	}
		    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    	}
	    }
	}else{ // already clean data
		string rawname=output_prefix;
		for(int j=0;j<infilelist1.size();j++){
			if(infilelist1.size() > 1)
				rawname = getfilename(infilelist1[j]);
			fprintf(stderr, "[ BatMeth2 ] Alignment %s, %s...\n", infilelist1[j].c_str(), infilelist2[j].c_str());
		    if(aligner == "BatMeth2")
		    	cmd="";
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + infilelist1[j] + " " + infilelist2[j] + " -t " + getstring(threads) + " --prefix " + outputdir + rawname + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + infilelist1[j] + " -b " + infilelist2[j] + " -d " + genome_index + " -o " + outputdir + rawname + ".sam -p " + getstring(threads); //+ " -r  0";
		    else if (aligner == "bismark2")
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters  + " -1 " + infilelist1[j] + " -2 " + infilelist2[j] + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + rawname + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    executeCMD(cmd.c_str(), outputdir, output_prefix);
		}
		if(aligner=="BatMeth2"){
            if(outformat=="BAM")
			    cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix1 + " -i " + input_prefix2 + " | samtools sort -m 2G -@ " + getstring(threads) + " -O BAM -o " + outputdir + output_prefix + ".sort.bam -";
            else
                cmd = abspath + "batmeth2-align" + " -L " + outputdir + output_prefix + ".align.log.txt" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix1 + " -i " + input_prefix2 + " -o " + outputdir + output_prefix + ".sam";
	    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
	    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    }else{
	    	if(infilelist1.size() > 1){
		    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
		    	for(int j=0; j< infilelist1.size(); j++){
		    		cmd = cmd + " " + outputdir + getfilename(infilelist1[j]) + ".sam";
		    	}
		    	cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
		    	executeCMD(cmd.c_str(), outputdir, output_prefix);
	    	}
	    }
	}

	
}

void annotation(string outputdir, string output_prefix){
    if(gfffile == "None")
        return;
    string methratio = outputdir + output_prefix + ".methratio.txt";
    string cmd;
    if(gfffile != "None")
        if(GTF)
            cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -gtf " + gfffile + " -m " + methratio + " -B -P --TSS --TTS --GENE -hs " + getstring(hs);
        else
            cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -gff " + gfffile + " -m " + methratio + " -B -P --TSS --TTS --GENE -hs " + getstring(hs);
    else if(bedfile != "None")
        cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -b " + bedfile + " -m " + methratio + " -B -P --TSS --TTS --GENE -hs "+ getstring(hs);
    else {
    	fprintf(stderr, "\nWarning: not defined gtf/gff/bed file, so skip annatation.\n");
    	return;
    }
    if(step != 0.01)
        cmd = cmd + " -s " + getstring(step);
    if(chromstep != 50000)
        cmd = cmd + " -S " + getstring(chromstep);
    if (coverage != 4)
        cmd = cmd + " -c " + getstring(coverage);
    if (maxcoverage != 1000)
        cmd = cmd + " -C " + getstring(maxcoverage);
    if (binCover != 1)
        cmd = cmd + " -nC " + getstring(binCover);
    if (distance != 2000)
        cmd = cmd + " -d " + getstring(distance);
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
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

//methylation data visulization
/*void visulize(){
    path1 = os.path.dirname(__file__)
    abspath = os.path.abspath(path1)
}*/

void printoutputfiles(string outputdir, string mkpath, string output_prefix){
	fprintf(stderr, "[ BatMeth2 ] Output files ...\n");
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

void mvpng(string outputdir, string mkpath, string output_prefix){
    vector<string> filelist;
    FileList(outputdir, output_prefix, "png", filelist);
    string srcstring = output_prefix + ".";
    for(int i = 0; i < filelist.size(); ++i){
        string oldpost = outputdir + "/" + filelist[i];
        string_replace(filelist[i], srcstring, "");
        string newpost = mkpath + "/images/" + filelist[i];
        fprintf(stderr, "move file %s to %s;\n", oldpost.c_str(), newpost.c_str());
        string cmd = "mv " + oldpost + " " + newpost;
        cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
}

void doc2html(string outputdir, string mkpath, string output_prefix){
	fprintf(stderr, "[ BatMeth2 ] Make batmeth2 report ...\n");
    string path = mkpath + "/images/";
    vector<string> filelist;
    FileList(path, "txt", "", filelist);
    for(int i = 0; i < filelist.size(); ++i){
        fprintf(stderr, "convert %s to html.\n", filelist[i].c_str());
        string cmd = "Rscript " + abspath + "doc2html.r " + outputdir + " " + output_prefix + " " + filelist[i];
        cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
}

void alignmentstate(string outputdir, string output_prefix, string mkpath){
    //string alignresults = outputdir + output_prefix + ".alignresults.txt";
    //FILE* Falignresults = File_Open(alignresults.c_str(), "w");
    //fprintf(Falignresults, "Value\tState\n");
    //fclose(Falignresults);
    if(outformat=="SAM"){
        string alignsortbam = "samtools sort -m 2G -@ " + getstring(threads) +" -o " + outputdir + output_prefix + ".sort.bam " + outputdir + output_prefix + ".sam";
        alignsortbam = alignsortbam + " >> " + outputdir + output_prefix + ".run.log 2>&1";
        executeCMD(alignsortbam.c_str(), outputdir, output_prefix);
    }
    //string alignsummarycmd = "samtools flagstat " + outputdir + output_prefix + ".sort.bam | tail -14 | sed 's/ /_/' | sed 's/ /_/' | sed 's/ /\\t/' | sed 's/_/ /g' > " + mkpath + output_prefix + ".alignresults.txt";
    //executeCMD(alignsummarycmd.c_str(), outputdir, output_prefix);
}

void methyPlot(string outputdir, string output_prefix ){
    string cmd="";
    //cmd = abspath + "GeneMethHeatmap " + outputdir + output_prefix + " None " + getstring(CG) + " " + getstring(CHG) + " " + getstring(CHH);
    //cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    //executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = abspath + "report2html " + "-p " + outputdir + output_prefix;
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    return;
    //
    cmd = abspath + "methyPlot " + outputdir + output_prefix + ".methBins.txt " + outputdir + output_prefix + ".Methygenome.pdf " + getstring(step) + " " + outputdir +output_prefix + ".AverMethylevel.txt " + outputdir + output_prefix + ".function.pdf TSS TTS " + outputdir +  output_prefix + ".Methylevel.txt " + outputdir +  output_prefix + ".Methenrich.pdf";
    //cmd = "{output_prefix}.bins..."
    //cmd = cmd.format(**locals())
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = "Rscript "+ abspath + "density_plot_with_methyl_oneSample_oneGff.r "+ outputdir + output_prefix + ".methBins.txt " + outputdir +  output_prefix + ".annoDensity.txt " + outputdir +  output_prefix + ".density.pdf " + output_prefix;
    //cmd = cmd.format(**locals())
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = "Rscript " + abspath + "mCdensity.r " + outputdir + output_prefix + ".mCdensity.txt " + outputdir + output_prefix + ".mCdensity.pdf " + outputdir + output_prefix + ".mCcatero.txt " + outputdir + output_prefix + ".mCcatero.pdf";
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = abspath + "GeneMethHeatmap " + outputdir + output_prefix + " None " + getstring(CG) + " " + getstring(CHG) + " " + getstring(CHH);
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
}

void visul2sample(){
    //string cmd = abspath + "GeneMethHeatmap {sample1_prefix} {sample2_prefix} {CG} {CHG} {CHG}";
}

void fastptrim(string outputdir, string output_prefix, string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean, bool pairedend){
    string prefix = outputdir + output_prefix;
    if(fastp != ""){
    	string cmd;
        if(pairedend)
            cmd = fastp + " -Y 0 -i " + input_prefix1 + " -I " + input_prefix2 + " -o " + input_clean1 + " -O " + input_clean2 + " -h " + prefix + ".html" \
                + " -j " + prefix + ".json";
        else
        	cmd = fastp + " -Y 0 -i " + input_prefix + " -o " + input_clean + " -h " + prefix + ".html" \
                + " -j " + prefix + ".json";
        cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
}

void mr2bw(string output_prefix, string outputdir){
    string cmd;
    cmd = "samtools faidx " + genome_index;
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = "python " + abspath + "batmeth2_to_bigwig.py -sort " + genome_index + ".fai " + output_prefix + ".methratio.txt";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    cmd = "python " + abspath + "batmeth2_to_bigwig.py -sort -strand " + genome_index + ".fai " + output_prefix + ".methratio.txt";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
}

void QCSingle(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix,
    string &clean_input){
    if(aligner == "no")
        return;
    if( (input_prefix == "None") || (output_prefix == "None")){
        fprintf(stderr, "Please check the pramater.\ngenome: %s\ninput: %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix.c_str(), output_prefix.c_str());
        exit(0);
    }

    //single-end QC
    std::vector <string> infilelist;
    std::vector<string> cleanfilelist;
    std::vector<string> outfilelist1;

    SplitString(input_prefix, infilelist, ",");

    string cmd="";
    //clean
    
    string input_clean;
    string fileformat;
    char temp[MAX_PATH];
    int i;
    string clenfiles="";
    string outfiles=""; //for trans fq files
    string cleanname = output_prefix;
    
    string out1_c2t;
    for(int j=0;j<infilelist.size();j++){
        for(i=0;i<infilelist[j].length();i++)
        temp[i] = infilelist[j][i];
        temp[i]='\0';
        get_fileformat(temp, fileformat);
        if( fileformat == "gz" || fileformat == "gzip" ){
            //gzfilelist[j] = 1;
            if(!cleanreads){
                input_clean = getfilename(string(temp)) + "clean.gz";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            //gzfilelist[j] = 0;
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
            if(clean_input == ""){
                clean_input = outputdir + cleanname;
            }else{
                clean_input = clean_input + "," + outputdir + cleanname;
            }
            infilelist[j]=cleanname;
            cleanfilelist.push_back(cleanname);
        }
        outfilelist1.push_back(out1_c2t);
    }

}

void QCPaired(string outputdir, string input_prefix, string input_prefix1, string input_prefix2, string output_prefix, 
    string &clean_input1, string &clean_input2){
    if(aligner == "no")
        return;
    if( (input_prefix1 == "") || (output_prefix == "") || (input_prefix2 == "")){
        fprintf(stderr, "\nError! Please check the pramater.\ngenome: %s\ninput: %s, %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix1.c_str(), input_prefix2.c_str(),output_prefix.c_str());
        exit(0);
    }

    //paired-end QC
    std::vector<string> cleanfilelist1;
    std::vector<string> cleanfilelist2;
    std::vector<string> outfilelist1;
    std::vector<string> outfilelist2;
    std::vector <string> infilelist1;
    std::vector <string> infilelist2;

    SplitString(input_prefix1, infilelist1, ",");
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

    string cleanname = output_prefix;
    for(int j=0;j<infilelist1.size();j++){
        for(i=0;i<infilelist1[j].length();i++)
            temp[i] = infilelist1[j][i];
        temp[i]='\0';
        get_fileformat(temp, fileformat);
        if( fileformat == "gz" || fileformat == "gzip" ) {
            //gzfilelist_1[j] = 1;
            if(!cleanreads){
                input_clean1 = getfilename(string(temp)) + "clean.gz";
                out1_c2t = getfilename(string(temp)) + "clean.c2t.fq.gz";
            }else{
                out1_c2t = getfilename(string(temp)) + "c2t.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            //gzfilelist_1[j] = 0;
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
            //gzfilelist_2[j] = 1;
            if(!cleanreads){
                input_clean2 = getfilename(string(temp)) + "clean.gz";
                out2_g2a = getfilename(string(temp)) + "clean.g2a.fq.gz";
            }else{
                out2_g2a = getfilename(string(temp)) + "g2a.fq.gz";
            }
        }else if(fileformat == "fq" || fileformat == "fastq"){
            //gzfilelist_2[j] = 0;
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
            if(clean_input1 == ""){
                clean_input1 = outputdir + input_clean1;
                clean_input2 = outputdir + input_clean2;
            }else{
                clean_input1 = clean_input1 + "," + outputdir + input_clean1;
                clean_input2 = clean_input2 + "," + outputdir + input_clean2;
            }
            infilelist1[j]=input_clean1;
            infilelist2[j]=input_clean2;
            cleanfilelist1.push_back(input_clean1);
            cleanfilelist2.push_back(input_clean2);
        }
        outfilelist1.push_back(out1_c2t);
        outfilelist2.push_back(out2_g2a);
        //alignment
        
    }

}

void alignment(string input_prefix1, string input_prefix2, string input_prefix, string outputdir,
    string output_prefix, bool pairedend){
    string clean_input1 = "", clean_input2 = "", clean_input = "";
    if(fastp!=""){
        fprintf(stderr, "[ BatMeth2 ] Clean reads ...\n");
        
        // read qc and c2t / g2a
        if(input_prefix1!="" && input_prefix2!=""){
            fprintf(stderr, "Process paired-end reads!\n");
            QCPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, clean_input1, clean_input2);
            // cleanfilelist1 cleanfilelist2
        }
        if(input_prefix!="None"){
            fprintf(stderr, "Process single-end reads\n");
            QCSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, clean_input);
            // cleanfilelist
        }
    }

    fprintf(stderr, "[ BatMeth2 ] Alignment ...\n");
    string align_result = outputdir + output_prefix + ".sort.bam";
    string cmd = abspath + "memalign c2t";
    if(fastp!=""){
        if(pairedend){
            cmd = cmd + " -1 " + clean_input1 + " -2 " + clean_input2;
        }
        if(input_prefix!="None"){
            cmd = cmd + " -i " + clean_input;
        }
    }else{
        if(pairedend){
            cmd = cmd + " -1 " + input_prefix1 + " -2 " + input_prefix2;
            //alignmentPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
        if(input_prefix!="None"){
            cmd = cmd + " -i " + input_prefix;
            //alignmentSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
    }
    
    cmd = cmd + " -o " + output_prefix; // + " -O " + outputdir; //for log outfile
    cmd = cmd + " | " + abspath + "bwame mem -t " + getstring(threads) + " -C -p -Y " + genome_index + ".batmeth2.fa -  | " \
      + "samtools sort -@ "+ getstring(threads) + " -o " + align_result + " - ";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
}

// whole pipeline for DNA methylation analysis. Contains alignment, calute meth level, DNA methylation annatation
// on gff file or bed region, DNA methylation visulization. Differentail analysis use diffmeth function.

void runpipe(string outputdir, string output_prefix, string mkpath, string input_prefix, string input_prefix1, string input_prefix2, bool pairedend){
	fprintf(stderr, "[ BatMeth2 ] Genome: %s\n", genome_index.c_str());
	fprintf(stderr, "[ BatMeth2 ] Annotation, gtf: %s; bed: %s;\n", gfffile.c_str(), bedfile.c_str());
	fprintf(stderr, "[ BatMeth2 ] Input file:  %s, %s %s\n", input_prefix.c_str(), input_prefix1.c_str(), input_prefix2.c_str());
	fprintf(stderr, "[ BatMeth2 ] Outfile prefix: %s\n", output_prefix.c_str());

    printparamter1(mkpath, input_prefix, input_prefix1, input_prefix2, outputdir, pairedend, output_prefix);
    printparamter2(mkpath, output_prefix);
    
    string clean_input1 = "", clean_input2 = "", clean_input = "";
    if(fastp!=""){
        fprintf(stderr, "[ BatMeth2 ] Clean reads ...\n");
        
        // read qc and c2t / g2a
        if(input_prefix1!="" && input_prefix2!=""){
            fprintf(stderr, "Process paired-end reads!\n");
            QCPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, clean_input1, clean_input2);
            // cleanfilelist1 cleanfilelist2
        }
        if(input_prefix!="None"){
            fprintf(stderr, "Process single-end reads\n");
            QCSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix, clean_input);
            // cleanfilelist
        }
    }

    fprintf(stderr, "[ BatMeth2 ] Alignment ...\n");
    string align_result = outputdir + output_prefix + ".sort.bam";
    string cmd = abspath + "memalign c2t";
    if(fastp!=""){
        if(pairedend){
            cmd = cmd + " -1 " + clean_input1 + " -2 " + clean_input2;
        }
        if(input_prefix!="None"){
            cmd = cmd + " -i " + clean_input;
        }
    }else{
        if(pairedend){
            cmd = cmd + " -1 " + input_prefix1 + " -2 " + input_prefix2;
            //alignmentPaired(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
        if(input_prefix!="None"){
            cmd = cmd + " -i " + input_prefix;
            //alignmentSingle(outputdir, input_prefix, input_prefix1, input_prefix2, output_prefix);
        }
    }
    
    cmd = cmd + " -o " + output_prefix; // + " -O " + outputdir; //for log outfile
    cmd = cmd + " | " + abspath + "bwame mem -t " + getstring(threads) + " -C -p -Y " + genome_index + ".batmeth2.fa -  | " \
      + "samtools sort -@ "+ getstring(threads) + " -o " + align_result + " - ";
    executeCMD(cmd.c_str(), outputdir, output_prefix);

    fprintf(stderr, "[ BatMeth2 ] Alignment summary ...\n");
    fprintf(stderr, "[ BatMeth2 ] Sorting align file ...\n");
    fprintf(stderr, "[ BatMeth2 ] ");
    alignmentstate(outputdir, output_prefix,  mkpath);
    fprintf(stderr, "[ BatMeth2 ] Methylation Status ...\n");
    align_result = outputdir + output_prefix + ".sort.bam";
    calmeth(align_result, outputdir, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Annotation ...\n");
    annotation(outputdir, output_prefix);
    string methratioLogfile = outputdir + output_prefix + ".log.txt";
    string newlogfile = mkpath + output_prefix + ".methbasic.txt";
    cmd = "cp ";
    cmd += methratioLogfile; cmd+=" ";
    cmd += newlogfile;
    cmd = cmd + " >> " + outputdir + output_prefix + ".run.log 2>&1";
    executeCMD(cmd.c_str(), outputdir, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] mr to bigwig ...\n");
    if(printbigwig)
        mr2bw(output_prefix, outputdir);
    fprintf(stderr, "[ BatMeth2 ] Visulization ...\n");
    methyPlot(outputdir, output_prefix );
    //mvpng(outputdir, mkpath, output_prefix);
    //printoutputfiles(outputdir, mkpath, output_prefix);
    //doc2html(outputdir, mkpath, output_prefix);
    fprintf(stderr, "[ BatMeth2 ] Done!\nBatMeth2 is a naive tools, if you meet any problems, please let us know. We will fix it asap!\nE-mail: qiangwei.zhou2013@gmail.com\n");
}

//"pipel", "index", "align", "calmeth", "anno", "visul", "diffmeth", "visuldiff", DManno"

void detect_mode(string mode, int Nparas, char* paramaters[], string outputdir, string input_prefix, \
  string input_prefix1, string input_prefix2, string output_prefix, string mkpath, bool pairedend, \
  string redss, string genome_index, string genomeprefix){
    string command = mode;
    string para = "";
    string tmp;
    if(command != "pipel"){
        for(int i=2; i< Nparas; i++){
            para += " ";
        	tmp = paramaters[i];
            para +=  tmp;
        }
    }
    string cmd;
    if(command == "index" || command == "index_rrbs")
        build_index(para, command, outputdir, output_prefix, redss, genome_index, genomeprefix);
    else if(command == "pipel")
        runpipe(outputdir, output_prefix, mkpath, input_prefix, input_prefix1, input_prefix2, pairedend);
    else if (command == "align")
        alignment(input_prefix1, input_prefix2, input_prefix, outputdir, output_prefix, pairedend);
    else if (command == "calmeth"){
        cmd = abspath + "calmeth " + para;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
    else if (command == "methyGff"){
        cmd = abspath + "methyGff " + para;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
    else if (command == "methyPlot"){
        cmd = abspath + "methyPlot" + para;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
    else if (command == "batDMR"){
        cmd = abspath + "batDMR" + para;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
    else if (command == "DMCplot"){
        cmd = abspath + "DMCannotationPlot" + para;
        executeCMD(cmd.c_str(), outputdir, output_prefix);
    }
    else if (command == "visul2sample")
        visul2sample();
    else{
        fprintf(stderr, "\ncan not detect any command mode!\n");
        usage();
        exit(0);
    }
}

