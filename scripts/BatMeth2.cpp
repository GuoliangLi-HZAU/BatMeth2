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

using std::string;
using std::vector;
//alignment
int threads = 8;
float mismatches = 0.05;
string genome_index = "";
bool pairedend=false;
string output_prefix = "None";

string input_prefix = "None";
string input_prefix1 = "";
string input_prefix2 = "";

string bismark2paramaters="";
string bsmapparamaters="";
string bwamethparamaters="";
// fastp
string fastp="";

//calmeth
int Qual = 10;
int redup = 0;
int region = 1000;
int sammeth = 0;
//calmeth and methyGff
int coverage = 5;
int maxcoverage = 1000;
int binCover = 3;
int chromstep = 50000;

//methyGff
float step = 0.025;
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
string mkpath;
string workdir;
bool cleanreads=true;
string programname;

void usage(){
    fprintf(stderr, "\nBatMeth2 [mode] [paramaters]\n");
    fprintf(stderr, "mode:  build_index, pipel, align, calmeth, annoation, methyPlot, batDMR, visul2sample, mkreport, DMCplot\n\n");
    fprintf(stderr, "Example:\n  BatMeth2 pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -1 in1.fq.gz -2 in2.fq.gz --gff gene.gff -f 1 -o outprefix\n");
    fprintf(stderr, "Or single-end:\n  BatMeth2 pipel --fastp ~/location/to/fastp -g genome_indexed_by_batmeth2 -i in.fq.gz --gff gene.gff -f 1 -o outprefix\n\n");
//    ~/software_devp/batmeth2/bin/test pipel --aligner bwa-meth --go ~/practice/Genome/arabidopsis/arabidopsis_bwa-meth/TAIR10_chr_all.fa --fastp ~/software/fastp/fastp -g ~/practice/Genome/arabidopsis/arabidopsis_batmeth2_index/TAIR10_chr_all.fa -1 test1.fq -2 test2.fq --gff ~/practice/Genome/arabidopsis/TAIR10.gene.modify.gff -f 1 -o pipel.clean
    fprintf(stderr, "\n[build_index]\n");
    fprintf(stderr, "    Usage: BatMeth2 build_index genomefile. (wgbs data, must run this step first), or\n");
    fprintf(stderr, "    Usage: BatMeth2 build_index rrbs genomefile. (rrbs data)\n");

    fprintf(stderr, "\n[pipel (Contains: align, calmeth, annoation, methyPlot, mkreport)]\n");
    fprintf(stderr, "[fastp location]\n");
    fprintf(stderr, "    --fastp    fastp program location.\n");
    fprintf(stderr, "    ** If --fastp is not defined, the input file should be clean data.\n");
    fprintf(stderr, "[select aligner]\n");
    fprintf(stderr, "    --aligner    BatMeth2(default), bwa-meth(v1), bsmap, bismark2, no (exit out_prefix.sam file, no need align again)\n");
    fprintf(stderr, "    --bismark2para    bismark2 paramaters, contained by \"paramaters\"\n");
    fprintf(stderr, "    --bsmsppara    bsmap paramaters\n");
    fprintf(stderr, "    --bwamethpara    bwameth paramaters\n");
    fprintf(stderr, "[other aligners paramaters]\n");
    fprintf(stderr, "    --go    Name of the genome, contaion index build by aligner. (bwa-meth/bismark2)\n");
    fprintf(stderr, "[main paramaters]\n");
    fprintf(stderr, "    -o [outprefix]    Name of output file prefix\n");
    fprintf(stderr, "    -O [out folder]    Output of result file to specified folder, default output to current folder (./)\n");
    fprintf(stderr, "[alignment paramaters]\n");
    fprintf(stderr, "    -i    Name of input file, support .fq/.fastq and .gz/.gzip format. if paired-end. please use -1, -2\n");
    fprintf(stderr, "    -1    Name of input file left end, if single-end. please use -i\n");
    fprintf(stderr, "    -2    Name of input file right end\n");
    fprintf(stderr, "    -i/-1/-2 can be comma-separated lists (no whitespace), only supported in BatMeth2 aligner.\n");
    fprintf(stderr, "    -g    Name of the genome mapped against, make sure build index first.\n");
    fprintf(stderr, "    -p <interger>    Launch <integer> threads\n");
    fprintf(stderr, "[calmeth paramaters]\n");
    fprintf(stderr, "    --Qual      calculate the methratio while read QulityScore >= Q. default:10\n");
    fprintf(stderr, "    --redup     REMOVE_DUP, 0 or 1, default 0\n");
    fprintf(stderr, "    --region    Bins for DMR calculate , default 1000bp .\n");
    fprintf(stderr, "    -f           for sam format outfile contain methState. [0 or 1], default: 0 (dont output this file).\n");
	fprintf(stderr, "    -n    maximum mismatches allowed due to seq. default 0.05 percentage of the read length. [0-0.3]\n");
    fprintf(stderr, "[calmeth and annoation paramaters]\n");
    fprintf(stderr, "    --coverage    >= <INT> coverage. default:5\n");
    fprintf(stderr, "    --binCover    >= <INT> nCs per region. default:3\n");
    fprintf(stderr, "    --chromstep   Chromsome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp)\n");
    fprintf(stderr, "[annoation paramaters]\n");
    fprintf(stderr, "    --gtf/--gff/--bed    Gtf or gff file / bed file\n");
    fprintf(stderr, "    --distance    DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000\n");
    fprintf(stderr, "    --step    Gene body and their flanking sequences using an overlapping sliding window of 0.05 of the sequence length at a step of 0.025 of the sequence length. So default step: 0.025 (2.5%)\n");
    fprintf(stderr, "    -C         <= <INT> coverage. default:1000\n");
    fprintf(stderr, "[MethyPlot paramaters]\n");
    fprintf(stderr, "    --CG       CG ratio for heatmap, [0-1], default 0.6\n");
    fprintf(stderr, "    --CHG      CHG ratio for heatmap, [0-1], default 0.2\n");
    fprintf(stderr, "    --CHH      CHH ratio for heatmap, [0-1], default 0.1\n");
    fprintf(stderr, "[mkreport paramaters]\n");
    fprintf(stderr, "    -o [outprefix]\n");
    fprintf(stderr, "    make a html report.\n");

    fprintf(stderr, "\n[align paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 align'\n");
    fprintf(stderr, "\n[calmeth paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 calmeth'\n");
    fprintf(stderr, "\n[annotion paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 annoation'\n");
    fprintf(stderr, "\n[methyPlot paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 methyPlot'\n");
    fprintf(stderr, "\n[batDMR paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 batDMR'\n");
    fprintf(stderr, "\n[visul2sample paramaters:]\n");
    fprintf(stderr, "    see the details in 'BatMeth2 visul2sample'\n\n");

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

void executeCMD(const char *cmd)
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

void printparamter1(string mkpath, string input_prefix, string input_prefix1, string input_prefix2, string outputdir){
	fprintf(stderr, "[ BatMeth2 ] Process Paramater file.\n");
    //paramater file
    string fparamater = mkpath + "/images/Paramater.txt";
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

void printparamter2(string mkpath){
	fprintf(stderr, "[ BatMeth2 ] Process Paramater file2.\n");
    //paramater2
    string fparamater2 = mkpath + "/images/Paramater2.txt";
    FILE* Fparamater = File_Open(fparamater2.c_str(), "w");
    fprintf(Fparamater, "Paramater\tValue\nThreads\t%d\nCalmeth\tparamaters\nQuality_Score\t%d\nredup\t%d\nmeth region length\t%d\nPrint methstate samfile\t%d\ncalmeth and methyGff\tparamaters\nCoverage\t%d\nmaxCoverage\t%d\nbinCoverage\t%d\nchromStep\t%d\nmethyGff\tParamaters\nGene bins step\t%.3f\nDistance of upstream and downstream\t%d  ", threads, Qual, redup, region, sammeth, coverage, maxcoverage, binCover, chromstep, step, distance);
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
void calmeth(string inputf);
void build_index(string para);
string get_path(string filepath);
void alignmentSingle();
void alignmentPaired();
void annoation();
void GetFileList(string PATH, FILE* outfile, string contain1);
void printoutputfiles();
void FileList(string PATH, string contain1, string contain2, vector<string>& files2);
void string_replace( std::string &strBig, const std::string &strsrc, const std::string &strdst);
void mvpng();
void doc2html();
void alignmentstate();
void mkreport();
void methyPlot();
void fastptrim(string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean);
void runpipe();
void detect_mode(string mode, int Nparas, char* paramaters[]);
void get_fileformat(char* processdir, string& processname);
int MAX_PATH = 1000;

std::vector <string> infilelist;
std::vector <string> infilelist1;
std::vector <string> infilelist2;
string outputdir="./";
int main(int argc, char* argv[])
{
	
	for(int i=1;i<argc;i++)
    {
        if(!strcmp(argv[i], "-i")){
        	input_prefix= argv[++i];
		    SplitString(input_prefix, infilelist, ",");
        }else if(!strcmp(argv[i], "-1")){
        	input_prefix1= argv[++i];
		    SplitString(input_prefix1, infilelist1, ",");
        }else if(!strcmp(argv[i], "-2")){
        	input_prefix2= argv[++i];
		    SplitString(input_prefix2, infilelist2, ",");
            pairedend=true;
        }else if(!strcmp(argv[i], "-o"))
            output_prefix= argv[++i];
        else if(!strcmp(argv[i], "-O")){
            outputdir= argv[++i];
            if(outputdir[outputdir.length()-1] != '/')
            	outputdir+="/";
            string cmd = "mkdir -p " + outputdir;
            executeCMD(cmd.c_str());
        }
        else if(!strcmp(argv[i], "-g"))
        	genome_index= argv[++i];
        else if(!strcmp(argv[i], "-p"))
            threads = atoi(argv[++i]);
        else if(!strcmp(argv[i], "-n"))
            mismatches = atof(argv[++i]);
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
        else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help") || !strcmp(argv[i], "-help")){
            usage();
            exit(0);
        }
    }

	if (argc < 2){
	    usage();
	    exit(0);
	}

    mode = argv[1];
    if(mode != "build_index" && mode != "pipel" && mode != "align" && mode != "calmeth" && mode != "annoation" && mode != "methyPlot" && mode != "batDMR" && 
    mode != "visul2sample" && mode != "mkreport" && mode != "DMCplot"){
    	fprintf(stderr, "\nNot a valid mode\n");
    	usage();
    	exit(0);
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
	fprintf(stderr, "[ Program directory ] %s\n[ Program name ] %s\n",abspathtmp, processname);
	programname = processname;

    if(mode == "align" && aligner == "BatMeth2"){
        if(argc < 4){ 
            string cmd = abspath + "batmeth2-align";
            executeCMD(cmd.c_str());
            exit(0);
        }
    }
    //for ann and bin file
    if(mode == "pipel" && genome_index == "" && genome_others != ""){
    	if(aligner == "BatMeth2"){
    		fprintf(stderr, "\nPlease defined genome_index location, when use batmeth2 aligner\n");
    		exit(0);
    	}
    	string cmd = abspath + "preGenome " + genome_others;
    	executeCMD(cmd.c_str());
    	genome_index = genome_others;
	}
    mkpath= outputdir + "/batmeth2_report_" + output_prefix;

    fprintf(stderr, "[ Workdir ] %s\n", workdirtmp);
    fprintf(stderr, "[ outputdir ] %s\n", outputdir.c_str());
	if (argc < 4){
	    if (mode == "pipel")
	        usage();
	    else
	        detect_mode(mode, argc, argv);
	    exit(0);
	}

	detect_mode(mode, argc, argv);

}

void calmeth(string inputf){
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
    if(coverage != 5)
        cmd = cmd + " -c " + getstring(coverage);
    if(binCover != 3)
        cmd = cmd + " -nC " + getstring(binCover);
    if(Qual != 10)
        cmd = cmd + " -Q " + getstring(Qual);
    if(region != 1000)
        cmd = cmd + " -R " + getstring(region);
    executeCMD(cmd.c_str());
    return;
}

void build_index(string para){
    if(para==""){
        fprintf(stderr, "Must have genome file.\n");
        exit(0);
    }
    string cmd = abspath + "build_all " + para; //sys.argv[2]
    executeCMD(cmd.c_str());
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

string get_path(string filepath)
{
	int n=filepath.find_last_of('/');
	string dirname=filepath.substr(0,n);
	//printf("\nss %s\n", dirname.c_str());
	return dirname;
}

void alignmentSingle(){
    if(aligner == "no")
        return;
    if((genome_index == "") || (input_prefix == "None") || (output_prefix == "None")){
        fprintf(stderr, "Please check the pramater.\ngenome: %s\ninput: %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix.c_str(), output_prefix.c_str());
        if(aligner == "BatMeth2"){
        	string cmd = abspath + programname;
            executeCMD(cmd.c_str());
        }
        exit(0);
    }
    string cmd="";
	if(!cleanreads){
		std::vector<string> cleanfilelist;
    	string input_clean;
    	string fileformat;
    	char temp[MAX_PATH];
    	int i;
    	string clenfiles="";
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
	    	string input_clean1;
	    	string input_clean2;
	    	fprintf(stderr, "[ BatMeth2 ] raw reads: %s; clean reads: %s\n", infilelist[j].c_str(), input_clean.c_str());
		    fastptrim(input_prefix1, input_prefix2, infilelist[j], input_clean1, input_clean2, outputdir + input_clean);
			infilelist[j]=input_clean;
			cleanfilelist.push_back(input_clean);
			//alignment
			fprintf(stderr, "[ BatMeth2 ] Alignment %s ...\n", input_clean.c_str());
		    if(aligner == "BatMeth2")
		    	if(clenfiles == "")
		    		clenfiles = outputdir + input_clean;
		    	else clenfiles = clenfiles + "," + outputdir + input_clean;
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + outputdir + input_clean + " -t " + getstring(threads) + " --prefix " + outputdir + input_clean + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + outputdir + input_clean + " -d " + genome_index + " -o " + outputdir + input_clean + ".sam -p " + getstring(threads);// + " -r 0";
		    else if(aligner == "bismark2") //need test get_path
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters + " -U " + outputdir + input_clean + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + input_clean + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    executeCMD(cmd.c_str());
	    }
	    if(aligner=="BatMeth2"){
	    	cmd = abspath + "batmeth2-align" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles + " -o " + outputdir + output_prefix + ".sam";
	    	executeCMD(cmd.c_str());
	    }else{
	    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
	    	for(int j=0; j< cleanfilelist.size(); j++){
	    		cmd = cmd + " " + outputdir + cleanfilelist[j] + ".sam";
	    	}
	    	executeCMD(cmd.c_str());
	    }
	}else{ // already clean data
		for(int j=0;j<infilelist.size();j++){
			fprintf(stderr, "[ BatMeth2 ] Alignment %s ...\n", infilelist[j].c_str());
		    if(aligner == "BatMeth2")
		    	cmd="";
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + infilelist[j] + " -t " + getstring(threads) + " --prefix " + outputdir + infilelist[j] + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + infilelist[j] + " -d " + genome_index + " -o " + outputdir + infilelist[j] + ".sam -p " + getstring(threads);// + " -r 0";
		    else if(aligner == "bismark2") //need test get_path
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters + " -U " + infilelist[j] + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + infilelist[j] + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    executeCMD(cmd.c_str());
		}
		if(aligner=="BatMeth2"){
	    	cmd = abspath + "batmeth2-align" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix + " -o " + outputdir + output_prefix + ".sam";
	    	executeCMD(cmd.c_str());
	    }else{
	    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
	    	for(int j=0; j< infilelist.size(); j++){
	    		cmd = cmd + " " + outputdir + infilelist[j] + ".sam";
	    	}
	    	executeCMD(cmd.c_str());
	    	string rmfile;
	    	//for(int j=0; j< infilelist.size(); j++){
	    	//	rmfile= infilelist[j] + ".sam";
	    	//	remove(rmfile)
	    	//}
	    }
	}


}

void alignmentPaired(){
    if(aligner == "no")
        return;
    if((genome_index == "") || (input_prefix1 == "") || (output_prefix == "None") || (input_prefix2 == "")){
        fprintf(stderr, "\nError! Please check the pramater.\ngenome: %s\ninput: %s, %s\noutput_prefix: %s\n", genome_index.c_str(), input_prefix1.c_str(), input_prefix2.c_str(),output_prefix.c_str());
        if(aligner == "BatMeth2"){
        	string cmd = abspath + programname;
            executeCMD(cmd.c_str());
        }
        exit(0);
    }
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
    	for(int j=0;j<infilelist1.size();j++){
	    	for(i=0;i<infilelist1[j].length();i++)
				temp[i] = infilelist1[j][i];
			temp[i]='\0';
	    	get_fileformat(temp, fileformat);
	    	if( fileformat == "gz" || fileformat == "gzip" ) 
	    		input_clean1 = string(temp) + "clean.gz";
	    	else if(fileformat == "fq" || fileformat == "fastq")
	    		input_clean1 = string(temp) + "clean.fq";
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
	    		input_clean2 = string(temp) + "clean.gz";
	    	else if(fileformat == "fq" || fileformat == "fastq")
	    		input_clean2 = string(temp) + "clean.fq";
	    	else{
	    		fprintf(stderr, "\n%s not a valid input files, should be fq/fastq or gz/gzip format!\n", fileformat.c_str());
	    		exit(0);
	    	}
	    	fprintf(stderr, "[ BatMeth2 ] raw reads: %s, %s; clean reads: %s, %s\n", infilelist1[j].c_str(), infilelist2[j].c_str(), input_clean1.c_str(), input_clean2.c_str());
		    fastptrim(infilelist1[j], infilelist2[j], input_prefix, outputdir + input_clean1, outputdir + input_clean2, input_clean);
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
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + outputdir + input_clean1 + " " + outputdir + input_clean2 + " -t " + getstring(threads) + " --prefix " + outputdir + input_clean1 + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + outputdir + input_clean1 + " -b " + outputdir + input_clean2 + " -d " + genome_index + " -o " + outputdir + input_clean1 + ".sam -p " + getstring(threads); //+ " -r  0";
		    else if (aligner == "bismark2")
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters  + " -1 " + outputdir + input_clean1 + " -2 " + outputdir + input_clean2 + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + input_clean1 + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    executeCMD(cmd.c_str());
		}
		if(aligner=="BatMeth2"){
		    cmd = abspath + "batmeth2-align" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + clenfiles1 + " -i " + clenfiles2 + " -o " + outputdir + output_prefix + ".sam";	    	
	    	executeCMD(cmd.c_str());
	    }else{
	    	cmd = "samtools merge -f -O SAM " + outputdir + output_prefix + ".sam";
	    	for(int j=0; j< cleanfilelist1.size(); j++){
	    		cmd = cmd + " " + outputdir + cleanfilelist1[j] + ".sam";
	    	}
	    	executeCMD(cmd.c_str());
	    }
	}else{ // already clean data
		for(int j=0;j<infilelist1.size();j++){
			fprintf(stderr, "[ BatMeth2 ] Alignment %s, %s...\n", infilelist1[j].c_str(), infilelist2[j].c_str());
		    if(aligner == "BatMeth2")
		    	cmd="";
		    else if(aligner == "bwa-meth")
		        cmd = "bwameth.py " + bwamethparamaters + " --reference " + genome_others + " " + infilelist1[j] + " " + infilelist2[j] + " -t " + getstring(threads) + " --prefix " + outputdir + infilelist1[j] + ".sam";
		    else if(aligner == "bsmap")
		        cmd = "bsmap" + bsmapparamaters + " -a " + infilelist1[j] + " -b " + infilelist2[j] + " -d " + genome_index + " -o " + outputdir + infilelist1[j] + ".sam -p " + getstring(threads); //+ " -r  0";
		    else if (aligner == "bismark2")
		        cmd = "bismark " + get_path(genome_others) + " " + bismark2paramaters  + " -1 " + infilelist1[j] + " -2 " + infilelist2[j] + " --bowtie2 " + " -p " + getstring(threads) + " --prefix " + outputdir + infilelist1[j] + " --sam";
		    else if (aligner == "no")
		        return;
		    else{
		        fprintf(stderr, "Please select correct aligner. (BatMeth2/bwa-meth/bsmap/bismark2)");
			    exit(0);
		    }
		    executeCMD(cmd.c_str());
		}
		if(aligner=="BatMeth2"){
			cmd = abspath + "batmeth2-align" + " -g " + genome_index + " -p " + getstring(threads) + " -i " + input_prefix1 + " -i " + input_prefix2 + " -o " + outputdir + output_prefix + ".sam";
	    	executeCMD(cmd.c_str());
	    }else{
	    	cmd = "samtools merge -f -O BAM " + outputdir + output_prefix + ".bam";
	    	for(int j=0; j< infilelist1.size(); j++){
	    		cmd = cmd + " " + outputdir + infilelist1[j] + ".sam";
	    	}
	    	executeCMD(cmd.c_str());
	    }
	}

	
}

void annoation(){
    if(gfffile == "None")
        return;
    string methratio = outputdir + output_prefix + ".methratio.txt";
    string cmd;
    if(gfffile != "None")
        if(GTF)
            cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -gtf " + gfffile + " -m " + methratio + " -B -P --TSS --TTS --GENE";
        else
            cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -gff " + gfffile + " -m " + methratio + " -B -P --TSS --TTS --GENE";
    else if(bedfile != "None")
        cmd = abspath + "methyGff" + " -o " + outputdir + output_prefix + " -G " + genome_index + " -b " + bedfile + " -m " + methratio + " -B -P --TSS --TTS --GENE";
    else {
    	fprintf(stderr, "\nWarning: not defined gtf/gff/bed file, so skip annatation.\n");
    	return;
    }
    if(step != 0.025)
        cmd = cmd + " -s " + getstring(step);
    if(chromstep != 50000)
        cmd = cmd + " -S " + getstring(chromstep);
    if (coverage != 5)
        cmd = cmd + " -c " + getstring(coverage);
    if (maxcoverage != 1000)
        cmd = cmd + " -C " + getstring(maxcoverage);
    if (binCover != 3)
        cmd = cmd + " -nC " + getstring(binCover);
    if (distance != 2000)
        cmd = cmd + " -d " + getstring(distance);
    executeCMD(cmd.c_str());
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

void printoutputfiles(){
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

void mvpng(){
    vector<string> filelist;
    FileList(outputdir, output_prefix, "png", filelist);
    string srcstring = output_prefix + ".";
    for(int i = 0; i < filelist.size(); ++i){
        string oldpost = outputdir + "/" + filelist[i];
        string_replace(filelist[i], srcstring, "");
        string newpost = mkpath + "/images/" + filelist[i];
        fprintf(stderr, "move file %s to %s;\n", oldpost.c_str(), newpost.c_str());
        string cmd = "mv " + oldpost + " " + newpost;
        executeCMD(cmd.c_str());
    }
}

void doc2html(){
	fprintf(stderr, "[ BatMeth2 ] Make batmeth2 report ...\n");
    string path = mkpath + "/images/";
    vector<string> filelist;
    FileList(path, "txt", "", filelist);
    for(int i = 0; i < filelist.size(); ++i){
        fprintf(stderr, "convert %s to html.\n", filelist[i].c_str());
        string cmd = "Rscript " + abspath + "doc2html.r " + outputdir + " " + output_prefix + " " + filelist[i];
        executeCMD(cmd.c_str());
    }
}

void alignmentstate(){
    string alignresults = outputdir + output_prefix + ".alignresults.txt";
    FILE* Falignresults = File_Open(alignresults.c_str(), "w");
    fprintf(Falignresults, "Value\tState\n");
    fclose(Falignresults);
    string alignsortbam = "samtools sort -@ " + getstring(threads) +" -o " + outputdir + output_prefix + ".sort.bam " + outputdir + output_prefix + ".sam";
    executeCMD(alignsortbam.c_str());
    string alignsummarycmd = "samtools flagstat " + outputdir + output_prefix + ".sort.bam >> " + alignresults + " && tail -14 " + alignresults + " | sed 's/ /_/' | sed 's/ /_/' | sed 's/ /\\t/' | sed 's/_/ /g' > " + mkpath + "/images/alignresults.txt";
    executeCMD(alignsummarycmd.c_str());
}

void mkreport(){
    if(output_prefix == "None"){
        fprintf(stderr, "\nplease define -o outfileprefix.\n");
        exit(0);
    }
    string copycmd="cp -r ";
    copycmd += abspath;
    copycmd += "../BatMeth2-Report " + outputdir + "batmeth2_report_";
    copycmd += output_prefix;
    executeCMD(copycmd.c_str());
    printparamter1(mkpath, input_prefix, input_prefix1, input_prefix2, outputdir);
    printparamter2(mkpath);
    alignmentstate();
    string methratioLogfile = output_prefix + output_prefix + ".log.txt";
    string newlogfile = mkpath + "/images/methratio.txt";
    string cmd = "cp ";
    cmd += methratioLogfile; cmd+=" ";
    cmd += newlogfile;
    executeCMD(cmd.c_str());
    mvpng();
    printoutputfiles();
    doc2html();
}

void methyPlot(){
    //
    string cmd = abspath + "methyPlot " + outputdir + output_prefix + ".methBins.txt " + outputdir + output_prefix + ".Methygenome.pdf " + getstring(step) + " " + outputdir +output_prefix + ".Methylevel.1.txt " + outputdir + output_prefix + ".function.pdf TSS TTS " + outputdir +  output_prefix + ".AverMethylevel.1.txt " + outputdir +  output_prefix + ".Methenrich.pdf";
    //cmd = "{output_prefix}.bins..."
    //cmd = cmd.format(**locals())
    executeCMD(cmd.c_str());
    cmd = "Rscript "+ abspath + "density_plot_with_methyl_oneSample_oneGff.r "+ outputdir + output_prefix + ".methBins.txt " + outputdir +  output_prefix + ".annoDensity.1.txt " + outputdir +  output_prefix + ".density.pdf " + output_prefix;
    //cmd = cmd.format(**locals())
    executeCMD(cmd.c_str());
    cmd = "Rscript " + abspath + "mCdensity.r " + outputdir + output_prefix + ".mCdensity.txt " + outputdir + output_prefix + ".mCdensity.pdf " + outputdir + output_prefix + ".mCcatero.txt " + outputdir + output_prefix + ".mCcatero.pdf";
    executeCMD(cmd.c_str());
    cmd = abspath + "GeneMethHeatmap " + outputdir + output_prefix + " None " + getstring(CG) + " " + getstring(CHG) + " " + getstring(CHH);
    executeCMD(cmd.c_str());
}

void visul2sample(){
    //string cmd = abspath + "GeneMethHeatmap {sample1_prefix} {sample2_prefix} {CG} {CHG} {CHG}";
}

void fastptrim(string input_prefix1, string input_prefix2, string input_prefix, string input_clean1, string input_clean2, string input_clean){
    if(fastp != ""){
    	string cmd;
        if(pairedend)
            cmd = fastp + " -Y 0 -i " + input_prefix1 + " -I " + input_prefix2 + " -o " + input_clean1 + " -O " + input_clean2;
        else
        	cmd = fastp + " -Y 0 -i " + input_prefix + " -o " + input_clean;
        executeCMD(cmd.c_str());
    }
}

// whole pipeline for DNA methylation analysis. Contains alignment, calute meth level, DNA methylation annatation
// on gff file or bed region, DNA methylation visulization. Differentail analysis use diffmeth function.

void runpipe(){
	fprintf(stderr, "[ BatMeth2 ] Genome: %s\n", genome_index.c_str());
	fprintf(stderr, "[ BatMeth2 ] Annotation, gtf: %s; bed: %s;\n", gfffile.c_str(), bedfile.c_str());
	fprintf(stderr, "[ BatMeth2 ] Input file:  %s, %s %s\n", input_prefix.c_str(), input_prefix1.c_str(), input_prefix2.c_str());
	fprintf(stderr, "[ BatMeth2 ] Outfile prefix: %s\n", output_prefix.c_str());

    string copycmd="cp -r ";
    copycmd += abspath;
    copycmd += "../BatMeth2-Report "+ outputdir + "batmeth2_report_";
    copycmd += output_prefix;
    executeCMD(copycmd.c_str());

    printparamter1(mkpath, input_prefix, input_prefix1, input_prefix2, outputdir);
    printparamter2(mkpath);
    fprintf(stderr, "[ BatMeth2 ] Alignment ...\n");
    string align_result = outputdir + output_prefix + ".sam";
    if(pairedend){
        alignmentPaired();
    }
    else{
        alignmentSingle();
    }
    fprintf(stderr, "[ BatMeth2 ] Alignment summary ...\n");
    fprintf(stderr, "[ BatMeth2 ] Sorting align file ...\n");
    fprintf(stderr, "[ BatMeth2 ] ");
    alignmentstate();
    fprintf(stderr, "[ BatMeth2 ] Methylation Status ...\n");
    align_result = outputdir + output_prefix + ".sort.bam";
    calmeth(align_result);
    fprintf(stderr, "[ BatMeth2 ] Annotation ...\n");
    annoation();
    string methratioLogfile = outputdir + output_prefix + ".log.txt";
    string newlogfile = mkpath + "/images/methratio.txt";
    string cmd = "cp ";
    cmd += methratioLogfile; cmd+=" ";
    cmd += newlogfile;
    executeCMD(cmd.c_str());
    fprintf(stderr, "[ BatMeth2 ] Visulization ...\n");
    methyPlot();
    mvpng();
    printoutputfiles();
    doc2html();
    fprintf(stderr, "[ BatMeth2 ] Done!\nBatMeth2 is a naive tools, if you meet any problems, please let us know. We will fix it asap!\nE-mail: qiangwei.zhou2013@gmail.com\n");
}

//"pipel", "index", "align", "calmeth", "anno", "visul", "diffmeth", "visuldiff", DManno"

void detect_mode(string mode, int Nparas, char* paramaters[]){
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
    if(command == "build_index")
        build_index(para);
    else if(command == "pipel")
        runpipe();
    else if (command == "align")
        if (pairedend)
            alignmentPaired();
        else
            alignmentSingle();
    else if (command == "calmeth"){
        cmd = abspath + "calmeth " + para;
        executeCMD(cmd.c_str());
    }
    else if (command == "annoation"){
        cmd = abspath + "methyGff " + para;
        executeCMD(cmd.c_str());
    }
    else if (command == "methyPlot"){
        cmd = abspath + "methyPlot" + para;
        executeCMD(cmd.c_str());
    }
    else if (command == "batDMR"){
        cmd = abspath + "batDMR" + para;
        executeCMD(cmd.c_str());
    }
    else if (command == "DMCplot"){
        cmd = abspath + "DMCannotationPlot" + para;
        executeCMD(cmd.c_str());
    }
    else if (command == "visul2sample")
        visul2sample();
    else if (command == "mkreport")
        mkreport();
    else{
        fprintf(stderr, "\ncan not detect any command mode!\n");
        usage();
        exit(0);
    }
}

