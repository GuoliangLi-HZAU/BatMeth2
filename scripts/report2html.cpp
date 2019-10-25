#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <string.h>
#include <sstream>

//#include <chrono>
//#include <ctime>
#include <time.h>

#define BATBUF 2000
#define BT2_VER "2.1"

using namespace std;

string command;
string list2string(double* list, int size);
string list2string_int(int* list, int size);
string list2string_long(long* list, int size);
void plot_chrom(ofstream& ofs, int x_len, double* chrom_cg_p, double* chrom_cg_n,
    double* chrom_chg_p, double* chrom_chg_n, double* chrom_chh_p, double* chrom_chh_n, char* chrom);
void printHeader(ofstream& ofs);
void printFooter(ofstream& ofs);
void outputRow(ofstream& ofs, string key, long v);
void outputRow(ofstream& ofs, string key, string v);
FILE* File_Open(const char* File_Name,const char* Mode);
string my_to_string(long n);
int main(int argc, char* argv[]){
	const char* Help_String="Command Format :  calmeth [options] -g GENOME  -i/-b <SamfileSorted/BamfileSorted> -m <methratio outfile prefix>\n"
		"\nUsage:\n"
		"\t-g|--genome           Genome\n"
		"\t-i|--input            Sam format file, sorted by chrom.\n"
		"\t-b|--binput           Bam format file, sorted by chrom.\n"
		//"\t-p|--threads          the number of threads.\n"
		"\t-n|--Nmismatch [float]  Number of mismatches, default 0.06 percentage of read length. [0-1]\n"
		"\t-m|--methratio        [MethFileNamePrefix]  Predix of methratio output file\n"
		"\t-Q [int]              caculate the methratio while read QulityScore >= Q. default:10\n"
		"\t-c|--coverage         >= <INT> coverage. default:5\n"
		"\t-nC		         >= <INT> nCs per region. default:5\n"
		"\t-R |--Regions         Bins for DMR caculate , default 1kb .\n"
		"\t--binsfile            DNA methylation level distributions in chrosome, default output file: {methratioPrefix}.methBins.txt\n"
		"\t-s|--step             Chrosome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp)\n"
		"\t-r|--remove_dup       REMOVE_DUP, default:false\n"
		"\t-f|--sam [outfile]    f for sam format outfile contain methState. default: sam format.\n"
		"\t--sam-seq-beforeBS    Converting BS read to the genome sequences.\n"
		"\t-h|--help";

    std::string enrichfile="";
    std::string htmlFile="batmeth2.html";
    string distri_mrfile="";
    string heatmap_cg_file="";
    string heatmap_chg_file="";
    string heatmap_chh_file="";
    string mrfile="";
    string prefix="";
    string chrom_mrfile="";
    char s2t[BATBUF];
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i], "-p") )
		{
			prefix=argv[++i];
		}
    }
    htmlFile = prefix + "." + htmlFile;
    enrichfile = prefix + ".AverMethylevel.1.txt";
    distri_mrfile = prefix + ".Methylevel.1.txt";
    heatmap_cg_file = prefix + ".1.txt.sorted.cg";
    heatmap_chg_file = prefix + ".1.txt.sorted.chg";
    heatmap_chh_file = prefix + ".1.txt.sorted.chh";
    mrfile = prefix + ".mCcatero.txt";
    chrom_mrfile = prefix + ".methBins.txt";
    string parafile = "Paramater.txt";
    string parafile2 = "Paramater2.txt";
    //awk -v FS="\t" -v OFS="\t" '{gsub(/ /,"_",$2);print $2,$1}' alignresults.txt > alignresults2.txt
    string alignfile = prefix + ".alignresults.txt";
    string mrbasicfile = prefix + ".methbasic.txt";

    ofstream ofs;
    ofs.open(htmlFile.c_str(), ifstream::out);
    printHeader(ofs);

    char tp_str[100]; //temp
    const unsigned int total = 200000;
    const unsigned int chrlen = 10000;
    std::string json_str = "";
    float mr=0;
    char* context = new char[8];

//----------------------------------------------------------------------------------------
//report header
    ofs << endl;
    ofs << "<h1 style='text-align:left;'><a href='https://github.com/GuoLiangLi-HZAU/BatMeth2' target='_blank' style='color:#663355;text-decoration:none;'> BatMeth2 Report</a>"<<endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n";
    ofs << "<div id='summary'>\n";

    //General
    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "batmeth2 version:", string(BT2_VER)+ " (<a href='https://github.com/GuoLiangLi-HZAU/BatMeth2'>https://github.com/GuoLiangLi-HZAU/BatMeth2</a>)");

    FILE* ParaFile = File_Open(parafile.c_str() , "r");
    char* opt1 = new char[1000];
    char* opt2 = new char[1000];
    while(fgets(s2t, BATBUF, ParaFile)!=0){
        if(s2t == NULL) continue;
        //sscanf(s2t, "%s%s", opt1, opt2);
        char *subarr = strtok(s2t,"\t");
        strcpy(opt1, subarr);
        subarr = strtok(NULL,"\t");
        strcpy(opt2, subarr);
        outputRow(ofs, opt1, opt2);
    }

    ofs << "</table>\n";
    ofs << "</div>\n";
    
    fclose(ParaFile);
    //===end General

    //Paramaters
    ofs << "<div class='subsection_title' onclick=showOrHide('paramaters')>Paramaters</div>\n";
    ofs << "<div id='paramaters'>\n";
    ofs << "<table class='summary_table'>\n";
 
    FILE* ParaFile2 = File_Open(parafile2.c_str() , "r");
    while(fgets(s2t, BATBUF, ParaFile)!=0){
        if(s2t == NULL) continue;
        char *subarr = strtok(s2t,"\t");
        strcpy(opt1, subarr);
        subarr = strtok(NULL,"\t");
        strcpy(opt2, subarr);
        outputRow(ofs, opt1, opt2);
    }

    ofs << "</table>\n";
    ofs << "</div>\n";
    
    fclose(ParaFile2);
    //===end Paramaters

    //=------ header
    ofs << "</div>\n";
    ofs << "</div>\n";

//-----------------------------------------------------------------------------------------
//alignment
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('alignment')><a name='summary'>Alignment</a></div>\n";
    ofs << "<div id='alignment'>\n";

    //align
    ofs << "<table class='summary_table'>\n";
 
    FILE* AlignF = File_Open(alignfile.c_str() , "r");
    while(fgets(s2t, BATBUF, AlignF)!=0){
        if(s2t == NULL) continue;
        char *subarr = strtok(s2t,"\t");
        strcpy(opt1, subarr);
        subarr = strtok(NULL,"\t");
        strcpy(opt2, subarr);
        outputRow(ofs, opt2, opt1);
    }

    ofs << "</table>\n";
    
    fclose(AlignF);
    //===end align

    ofs << "</div>\n";

//-----------------------------------------------------------------------------------------
//mrbasic
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('mrbasic')><a name='summary'>DNA Methylation level</a></div>\n";
    ofs << "<div id='mrbasic'>\n";

    //align
    ofs << "<table class='summary_table'>\n";
 
    FILE* MrbasicF = File_Open(mrbasicfile.c_str() , "r");
    while(fgets(s2t, BATBUF, MrbasicF)!=0){
        if(s2t == NULL) continue;
        //sscanf(s2t, "%s%s", opt1, opt2);
        char *subarr = strtok(s2t,"\t");
        strcpy(opt1, subarr);
        subarr = strtok(NULL,"\t");
        strcpy(opt2, subarr);
        outputRow(ofs, opt1, opt2);
    }

    ofs << "</table>\n";
    
    delete []opt1;
    delete []opt2;
    fclose(MrbasicF);
    //===end align

    ofs << "</div>\n";

//-----------------------------------------------------------------------------------------
//DNA methylation level
    // cater_mr header
    ofs << "<div class='section_div'>\n";

    //html of function part
    ofs << "<div id='cater_mr_figure'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('cater_mr')><a name='summary'>DNA methylation level</a></div>\n";
    ofs << "<div id='cater_mr' style='display:flex;'>\n";
    ofs << "<div class='col-md-4 figure' id='plot_mc' style='height:400px;width:500px;display:flex;'></div>\n";
    ofs << "<div class='col-md-4 figure' id='plot_cater_mr' style='height:400px;width:680px;display:flex;'></div>\n";
    ofs << "</div>\n";

    // script of plotly
    ofs << "\n<script type=\"text/javascript\">" << endl;
    json_str = "";

    FILE* Cater_MR = File_Open(mrfile.c_str(), "r");

    long* allc = new long[5];
    long* cpg = new long[5];
    char* catero = new char[16];
    double* mC = new double[4];
    long count = 0;
    while(fgets(s2t, BATBUF, Cater_MR)!=0){
        if(s2t == NULL) continue;
        if(s2t[0] != 'm'){
            sscanf(s2t, "%s%ld", catero, &count);
            
            if(strcmp(catero, "M")==0){
                allc[0] = count;
            }else if(strcmp(catero, "Mh")==0){
                allc[1] = count;
            }else if(strcmp(catero, "H")==0){
                allc[2] = count;
            }else if(strcmp(catero, "hU")==0){
                allc[3] = count;
            }else if(strcmp(catero, "U")==0){
                allc[4] = count;
            }else if(strcmp(catero, "CpG_M")==0){
                cpg[0] = count;
            }else if(strcmp(catero, "CpG_Mh")==0){
                cpg[1] = count;
            }else if(strcmp(catero, "CpG_H")==0){
                cpg[2] = count;
            }else if(strcmp(catero, "CpG_hU")==0){
                cpg[3] = count;
            }else if(strcmp(catero, "CpG_U")==0){
                cpg[4] = count;
            }
        }else{
            sscanf(s2t, "%s%f", catero, &mr);
            
            if(strcmp(catero, "mC")==0){
                mC[0] = mr;
            }else if(strcmp(catero, "mCG")==0){
                mC[1] = mr;
            }else if(strcmp(catero, "mCHG")==0){
                mC[2] = mr;
            }else if(strcmp(catero, "mCHH")==0){
                mC[3] = mr;
            }
        }

    }
    delete []catero;
    fclose(Cater_MR);

    //-- meth bar
    json_str += "var trace_mc ={";
    json_str += "x: ['mC', 'mCG', 'mCHG', 'mCHH'],";
    json_str += "y: [" + list2string(mC, 4) + "],";
    //json_str += "name: 'allc',";
    json_str += "type:'bar',";
    json_str += "marker: {color: '#6495ED'}";
    json_str += "};\n";

    json_str += "var data = [trace_mc];\n";
    json_str += "var layout={yaxis: { title:'DNA methylation level'} };\n";

    json_str += "Plotly.newPlot('plot_mc', data, layout);\n";

    //-- catero bar
    json_str += "var trace_allc ={";
    json_str += "x: ['M (>0.8)', 'Mh (>0.6)', 'H (>0.4)', 'hU (>0.2)', 'U (<=0.2)'],";
    json_str += "y: [" + list2string_long(allc, 5) + "],";
    json_str += "name: 'allc',";
    json_str += "type:'bar',";
    json_str += "marker: {color: '#1C86EE'}";
    json_str += "};\n";

    json_str += "var trace_cpg ={";
    json_str += "x: ['M (>0.8)', 'Mh (>0.6)', 'H (>0.4)', 'hU (>0.2)', 'U (<=0.2)'],";
    json_str += "y: [" + list2string_long(cpg, 5) + "],";
    json_str += "name: 'cpg',";
    json_str += "type:'bar',";
    json_str += "marker: {color: '#EEAEEE'}";
    json_str += "};\n";

    json_str += "var data = [trace_allc, trace_cpg];\n";
    json_str += "var layout={yaxis: { title:'The number of Methylation sites'} };\n";

    json_str += "Plotly.newPlot('plot_cater_mr', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;

    //distri footer
    ofs << "</div>\n";
    ofs << "</div>\n";

    delete []allc;
    delete []cpg;
    delete []mC;

//-----------------------------------------------------------------------------------------
//DNA methylation level in chromosome
    char* chrom = new char[100];
    char* old_chr = new char[100];
    strcpy(old_chr, "nullchrom");
    int chr_loci=0; int old_chr_loci=0;
    FILE* Chrom_MR = File_Open(chrom_mrfile.c_str(), "r");
    // chrom_mr header
    ofs << "<div class='section_div'>\n";

    //html of function part
    ofs << "<div id='chrom_mr_figure'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('chrom_mr')><a name='summary'>DNA methylation level in chromosome</a></div>\n";
    int ploting_chr = 0; int plotchr=0; //plot number of chromosome
    string plot_name="";
    
    while(fgets(s2t, BATBUF, Chrom_MR)!=0){
        if(s2t == NULL) continue;
        sscanf(s2t, "%s%d%f%s", chrom, &chr_loci, &mr, context);
        if(strcmp(old_chr, chrom)!=0){
            if(strcmp(old_chr, "nullchrom")!=0 && old_chr_loci>20){
                plot_name=old_chr;
                if(ploting_chr <= plotchr){
                    ofs << "<div id='chrom_mr' style='display:flex;'>\n";
                    ofs << "<div class='col-md-4 figure' id='plot_chrom_cg_" + plot_name + "' style='display:flex;height:400px;width:600px;'></div>\n";
                    ofs << "<div class='col-md-4 figure' id='plot_chrom_chg_" + plot_name + "' style='display:flex;height:400px;width:600px;'></div>\n";
                    ofs << "<div class='col-md-4 figure' id='plot_chrom_chh_" + plot_name + "' style='display:flex;height:400px;width:600px;'></div>\n";
                    ofs << "</div>\n";
                }
                ploting_chr++;
            }
            strcpy(old_chr, chrom);
        }
        old_chr_loci = chr_loci;   
    }
    
    rewind(Chrom_MR);
    double* chrom_cg_p = new double[chrlen];
    double* chrom_chg_p = new double[chrlen];
    double* chrom_chh_p = new double[chrlen];
    double* chrom_cg_n = new double[chrlen];
    double* chrom_chg_n = new double[chrlen];
    double* chrom_chh_n = new double[chrlen];
    for(int i=0;i<chrlen; i++){
        chrom_cg_p[i] = 0;
        chrom_cg_n[i] = 0;
        chrom_chg_p[i] = 0;
        chrom_chg_n[i] = 0;
        chrom_chh_p[i] = 0;
        chrom_chh_n[i] = 0;
    }
    unsigned int chrom_x = 0;
    ploting_chr=0;
    strcpy(old_chr, "nullchrom");
    chr_loci=0; old_chr_loci=0;
    mr=0;
    while(fgets(s2t, BATBUF, Chrom_MR)!=0){
        if(s2t == NULL) continue;
        sscanf(s2t, "%s%d%f%s", chrom, &chr_loci, &mr, context);
        if(strcmp(context, "CG")==0){
            if(mr>0)  chrom_cg_p[chr_loci-1] = mr;
            else if(mr < 0) chrom_cg_n[chr_loci-1] = mr;
            else{
                chrom_cg_p[chr_loci-1] = mr;
                chrom_cg_n[chr_loci-1] = mr;
            }
        }else if(strcmp(context, "CHG")==0){
            if(mr>0)  chrom_chg_p[chr_loci-1] = mr;
            else if(mr < 0) chrom_chg_n[chr_loci-1] = mr;
            else{
                chrom_chg_p[chr_loci-1] = mr;
                chrom_chg_n[chr_loci-1] = mr;
            }
        }else if(strcmp(context, "CHH")==0){
            if(mr>0)  chrom_chh_p[chr_loci-1] = mr;
            else if(mr < 0) chrom_chh_n[chr_loci-1] = mr;
            else{
                chrom_chh_p[chr_loci-1] = mr;
                chrom_chh_n[chr_loci-1] = mr;
            }
        }

        if(strcmp(old_chr, chrom)!=0){
            if(strcmp(old_chr, "nullchrom")!=0 & old_chr_loci>20){
                //plot
                chrom_cg_p[old_chr_loci]='\0';
                chrom_cg_n[old_chr_loci]='\0';
                chrom_chg_p[old_chr_loci]='\0';
                chrom_chg_n[old_chr_loci]='\0';
                chrom_chh_p[old_chr_loci]='\0';
                chrom_chh_n[old_chr_loci]='\0';
                if(ploting_chr <= plotchr){
                    plot_chrom(ofs, old_chr_loci, chrom_cg_p, chrom_cg_n, chrom_chg_p, chrom_chg_n, chrom_chh_p, chrom_chh_n, old_chr);
                }
                ploting_chr++;
            }
            strcpy(old_chr, chrom);
        }
        old_chr_loci = chr_loci;
    }
    delete []old_chr;
    delete []chrom;
    delete []chrom_cg_p;
    delete []chrom_chg_p;
    delete []chrom_chh_p;
    delete []chrom_cg_n;
    delete []chrom_chg_n;
    delete []chrom_chh_n;
    fclose(Chrom_MR);

    //footer
    ofs << "</div>\n";
    ofs << "</div>\n";

//-------------------------------------------------------------------------------------------------------------------------
//DNA methylation distribution
    // distri_mr header
    ofs << "<div class='section_div'>\n";

    //html of function part
    ofs << "<div id='distri_mr_figure'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('distri_mr')><a name='summary'>DNA methylation distribution</a></div>\n";
    ofs << "<div id='distri_mr'>\n";
    ofs << "<div class='col-md-4 figure' id='plot_distri_mr' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    // script of plotly
    ofs << "\n<script type=\"text/javascript\">" << endl;

    FILE* Distri_MR = File_Open(distri_mrfile.c_str(), "r");
    int i=0;
    double** cg_dis_mr = new double*[3];
    for(int i=0;i<3;i++)
        cg_dis_mr[i] = new double[total];
    int j_len=0;
    while(fgets(s2t, BATBUF, Distri_MR)!=0){
        char *subarr = strtok(s2t,"\t");
        j_len=0;
        while(subarr!=NULL){
            if(subarr[0]=='C') {
                subarr = strtok(NULL,"\t");
                continue;
            }
            cg_dis_mr[i][j_len] = atof(subarr);
            
            subarr = strtok(NULL,"\t");
            j_len++;
        }
        i++;
    }
    fclose(Distri_MR);

    //html
    int* x_axis = new int[j_len];
    for(int i=0;i<j_len+1;i++){
        x_axis[i]=i+1;
    }
    json_str += "var trace_distri_cg ={";
    json_str += "x: [" + list2string_int(x_axis, j_len) +"],";
    json_str += "y: [" + list2string(cg_dis_mr[0], j_len) + "],";
    json_str += "name: 'cg',";
    json_str += "mode:'lines',";
    json_str += "marker: {color: '#1C86EE'}";
    json_str += "};\n";

    json_str += "var trace_distri_chg ={";
    json_str += "x: [" + list2string_int(x_axis, j_len) +"],";
    json_str += "y: [" + list2string(cg_dis_mr[1], j_len) + "],";
    json_str += "name: 'chg',";
    json_str += "mode:'lines',";
    json_str += "marker: {color: '#FF4136'}";
    json_str += "};\n";

    json_str += "var trace_distri_chh ={";
    json_str += "x: [" + list2string_int(x_axis, j_len) +"],";
    json_str += "y:[" + list2string(cg_dis_mr[2], j_len) + "],";
    json_str += "name: 'chh',";
    json_str += "mode:'lines',";
    json_str += "marker:{color:'#FF851B'}";
    json_str += "};\n";

    json_str += "var data = [trace_distri_cg, trace_distri_chg, trace_distri_chh];\n";
    json_str += "var layout={yaxis: { title:'DNA Methylation level(%)'}, xaxis: {title: 'Function elements'} };\n";

    json_str += "Plotly.newPlot('plot_distri_mr', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;

    //distri footer
    ofs << "</div>\n";
    ofs << "</div>\n";
    for(int i=0;i<3;i++)
        delete []cg_dis_mr[i];
    delete []cg_dis_mr;
//----------------------------------------------------------------------------------------------------------
//DNA methylation heatmap
    // header
    ofs << "<div class='section_div'>\n";

    //html of function part
    ofs << "<div id='heatmap_cg_figure'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('heatmap')><a name='summary'>Heatmap of DNA methylation in function elements</a></div>\n";
    ofs << "<div id='heatmap' style='display:flex'>\n";
    ofs << "<div class='col-md-4 figure' id='plot_heatmap_cg' style='display:flex;height:500px; width: 400px'></div>\n";
    ofs << "<div class='col-md-4 figure' id='plot_heatmap_chg' style='display:flex;height:500px; width: 400px'></div>\n";
    ofs << "<div class='col-md-4 figure' id='plot_heatmap_chh' style='display:flex;height:500px; width: 400px'></div>\n";
    ofs << "</div>\n";

    // script of plotly
    ofs << "\n<script type=\"text/javascript\">" << endl;
    json_str = "";
//Empty Value Containers
    json_str += "var xValues = [];";
    json_str += "var yValues = [];";
    json_str += "var yShift = [];";
    json_str += "var finalX = [];";
    json_str += "var finalY = [];\n";
    json_str += "var data = [\n{\n z: [";

    FILE* Heatmap_Cg = File_Open(heatmap_cg_file.c_str(), "r");
    double* cg_heatmap = new double[total];
    j_len=0; int x_ax=0; 
    while(fgets(s2t, BATBUF, Heatmap_Cg)!=0){
        char *subarr = strtok(s2t,"\t");
        j_len=0;
        while(subarr!=NULL){
            if(subarr[0]!='0') {
                subarr = strtok(NULL,"\t");
                continue;
            }
            //printf("%s\n",subarr);
            mr=atof(subarr);
            if(mr>0.6) mr=0.6;
            cg_heatmap[j_len] = mr;
            subarr = strtok(NULL,"\t");
            j_len++;
        }
        // ---------- html part
        int* x_axis;
        if(x_ax==0){
            x_axis = new int[j_len];
            for(int i=0;i<j_len+1;i++){
                x_axis[i]=i+1;
            }
        }
        if(x_ax==0)
            json_str += "[" + list2string(cg_heatmap, j_len) +"]";
        else
            json_str += ", [" + list2string(cg_heatmap, j_len) +"]";
        x_ax++;
    }
    json_str += "],\n";
    json_str += "colorscale : 'Portland',";
    json_str += "type: 'heatmap' }";
    json_str += "];";

    json_str += "var layout = {";
    json_str += "title: 'CG',";
    json_str += "xaxis: {visible: false},";
    json_str += "yaxis: {visible: false}";
    json_str += "};\n";

    fclose(Heatmap_Cg);
    delete []cg_heatmap;

    json_str += "Plotly.newPlot('plot_heatmap_cg', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    //==================================================================
    ofs << "\n<script type=\"text/javascript\">" << endl;
    json_str = "";
//Empty Value Containers
    json_str += "var xValues = [];";
    json_str += "var yValues = [];";
    json_str += "var yShift = [];";
    json_str += "var finalX = [];";
    json_str += "var finalY = [];\n";
    json_str += "var data = [\n{\n z: [";

    FILE* Heatmap_Chg = File_Open(heatmap_chg_file.c_str(), "r");
    double* chg_heatmap = new double[total];
    j_len=0; x_ax=0; 
    while(fgets(s2t, BATBUF, Heatmap_Chg)!=0){
        char *subarr = strtok(s2t,"\t");
        j_len=0;
        while(subarr!=NULL){
            if(subarr[0]!='0') {
                subarr = strtok(NULL,"\t");
                continue;
            }
            //printf("%s\n",subarr);
            mr=atof(subarr);
            if(mr>0.1) mr=0.1;
            chg_heatmap[j_len] = mr;
            subarr = strtok(NULL,"\t");
            j_len++;
        }
        // ---------- html part
        int* x_axis;
        if(x_ax==0){
            x_axis = new int[j_len];
            for(int i=0;i<j_len+1;i++){
                x_axis[i]=i+1;
            }
        }
        if(x_ax==0)
            json_str += "[" + list2string(chg_heatmap, j_len) +"]";
        else
            json_str += ", [" + list2string(chg_heatmap, j_len) +"]";
        x_ax++;
    }
    json_str += "],\n";
    json_str += "colorscale : 'Portland',";
    json_str += "type: 'heatmap' }";
    json_str += "];";

    json_str += "var layout = {";
    json_str += "title: 'CHG',";
    json_str += "xaxis: {visible: false},";
    json_str += "yaxis: {visible: false}";
    json_str += "};\n";

    fclose(Heatmap_Chg);
    delete []chg_heatmap;

    json_str += "Plotly.newPlot('plot_heatmap_chg', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    //====================================================================
    ofs << "\n<script type=\"text/javascript\">" << endl;
    json_str = "";
//Empty Value Containers
    json_str += "var xValues = [];";
    json_str += "var yValues = [];";
    json_str += "var yShift = [];";
    json_str += "var finalX = [];";
    json_str += "var finalY = [];\n";
    json_str += "var data = [\n{\n z: [";

    FILE* Heatmap_Chh = File_Open(heatmap_chh_file.c_str(), "r");
    double* chh_heatmap = new double[total];
    j_len=0; x_ax=0; 
    while(fgets(s2t, BATBUF, Heatmap_Chh)!=0){
        char *subarr = strtok(s2t,"\t");
        j_len=0;
        while(subarr!=NULL){
            if(subarr[0]!='0') {
                subarr = strtok(NULL,"\t");
                continue;
            }
            //printf("%s\n",subarr);
            mr=atof(subarr);
            if(mr>0.1) mr=0.1;
            chh_heatmap[j_len] = mr;
            subarr = strtok(NULL,"\t");
            j_len++;
        }
        // ---------- html part
        int* x_axis;
        if(x_ax==0){
            x_axis = new int[j_len];
            for(int i=0;i<j_len+1;i++){
                x_axis[i]=i+1;
            }
        }
        if(x_ax==0)
            json_str += "[" + list2string(chh_heatmap, j_len) +"]";
        else
            json_str += ", [" + list2string(chh_heatmap, j_len) +"]";
        x_ax++;
    }
    json_str += "],\n";
    json_str += "colorscale : 'Portland',";
    json_str += "type: 'heatmap' }";
    json_str += "];";

    json_str += "var layout = {";
    json_str += "title: 'CHH',";
    json_str += "xaxis: {visible: false},";
    json_str += "yaxis: {visible: false}";
    json_str += "};\n";

    fclose(Heatmap_Chh);
    delete []chh_heatmap;

    json_str += "Plotly.newPlot('plot_heatmap_chh', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    //====================================================================
    ofs << "</div>\n";
    ofs << "</div>\n";
    
//----------------------------------------------------------------------------------------------------------
//Enrichment of DNA methylation
    // enrich header
    ofs << "<div class='section_div'>\n";
    
    //html of function part
    ofs << "<div id='enrich_figure'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('enrichplot')><a name='summary'>Enrichment of DNA methylation</a></div>\n";
    ofs << "<div id='enrichplot' style='display:flex'>\n";
    ofs << "<div class='col-md-4 figure' id='plot_enrich_cg' style='height:500px;width:400px;'></div>\n";
    ofs << "<div class='col-md-4 figure' id='plot_enrich_chg' style='height:500px;width:400px;'></div>\n";
    ofs << "<div class='col-md-4 figure' id='plot_enrich_chh' style='height:500px;width:400px;'></div>\n";
    ofs << "</div>\n";

    // script of plotly
    ofs << "\n<script type=\"text/javascript\">" << endl;

    json_str = "";

    double* cg_down_mr = new double[total];
    double* cg_body_mr = new double[total];
    double* cg_up_mr = new double[total];
    double* chg_down_mr = new double[total];
    double* chg_body_mr = new double[total];
    double* chg_up_mr = new double[total];
    double* chh_down_mr = new double[total];
    double* chh_body_mr = new double[total];
    double* chh_up_mr = new double[total];
    unsigned int cg_1=0;
    unsigned int cg_2=0;
    unsigned int cg_3=0;
    unsigned int chg_1=0;
    unsigned int chg_2=0;
    unsigned int chg_3=0;
    unsigned int chh_1=0;
    unsigned int chh_2=0;
    unsigned int chh_3=0;
    //temp
    char* region = new char[8];
    mr=0;
    FILE* Enrich_File=File_Open(enrichfile.c_str(),"r");
    while( fgets(s2t, BATBUF, Enrich_File)!=0 ){
        if(s2t!=NULL && s2t[1] == 'o') continue;
        sscanf(s2t, "%s%s%f", context, region, &mr);

        sprintf( tp_str, "%0.2f", mr); //截取后转化为字符串类型
        sscanf( tp_str, "%f", &mr);
        if(strcmp(context, "CG")==0){
            if(strcmp(region, "UP")==0){
                cg_up_mr[cg_1] = mr;
                cg_1++;
            }else if(strcmp(region, "BODY")==0){
                cg_body_mr[cg_2] = mr;
                cg_2++;
            }else if(strcmp(region, "DOWN")==0){
                cg_down_mr[cg_3] = mr;
                cg_3++;
            }
        }else if(strcmp(context, "CHG")==0){
            if(strcmp(region, "UP")==0){
                chg_up_mr[chg_1] = mr;
                chg_1++;
            }else if(strcmp(region, "BODY")==0){
                chg_body_mr[chg_2] = mr;
                chg_2++;
            }else if(strcmp(region, "DOWN")==0){
                chg_down_mr[chg_3] = mr;
                chg_3++;
            }
        }else if(strcmp(context, "CHH")==0){
            if(strcmp(region, "UP")==0){
                chh_up_mr[chh_1] = mr;
                chh_1++;
            }else if(strcmp(region, "BODY")==0){
                chh_body_mr[chh_2] = mr;
                chh_2++;
            }else if(strcmp(region, "DOWN")==0){
                chh_down_mr[chh_3] = mr;
                chh_3++;
            }
        }
    }
    fclose(Enrich_File);
    //html
    json_str += "var trace_cg_up ={";
    json_str += "y: [" + list2string(cg_up_mr, cg_1) + "],";
    json_str += "name: 'cg_upper',";
    json_str += "type:'box',";
    json_str += "marker: {color: '#1C86EE'}";
    json_str += "};\n";

    json_str += "var trace_cg_body ={";
    json_str += "y: [" + list2string(cg_body_mr, cg_2) + "],";
    json_str += "name: 'cg_body',";
    json_str += "type:'box',";
    json_str += "marker: {color: '#FF4136'}";
    json_str += "};\n";

    json_str += "var trace_cg_down ={";
    json_str += "y:[" + list2string(cg_down_mr, cg_1) + "],";
    json_str += "name: 'cg_down',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#FF851B'}";
    json_str += "};\n";

    json_str += "var data = [trace_cg_up, trace_cg_body, trace_cg_down];\n";
    json_str += "var layout={title: 'CG', xaxis: {tickangle: 45}, yaxis: { title:'Methylation level in CG (%)', zeroline: false}, boxmode: 'group' };\n";

    json_str += "Plotly.newPlot('plot_enrich_cg', data, layout);\n";
    // chg
    json_str += "var trace_chg_up ={";
    json_str += "y:[" + list2string(chg_up_mr, chg_1) + "],";
    json_str += "name: 'chg_upper',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#1C86EE'}";
    json_str += "};\n";

    json_str += "var trace_chg_body ={";
    json_str += "y:[" + list2string(chg_body_mr, chg_2) + "],";
    json_str += "name: 'chg_body',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#FF4136'}";
    json_str += "};\n";

    json_str += "var trace_chg_down ={";
    json_str += "y:[" + list2string(chg_down_mr, chg_1) + "],";
    json_str += "name: 'chg_down',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#FF851B'}";
    json_str += "};\n";

    json_str += "var data = [trace_chg_up, trace_chg_body, trace_chg_down];\n";
    json_str += "var layout={title: 'CHG', xaxis: {tickangle: 45}, yaxis: { title:'Methylation level in CHG (%)', zeroline: false}, boxmode: 'group' };\n";

    json_str += "Plotly.newPlot('plot_enrich_chg', data, layout);\n";
    // chh
    json_str += "var trace_chh_up ={";
    json_str += "y:[" + list2string(chh_up_mr, chh_1) + "],";
    json_str += "name: 'chh_upper',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#1C86EE'}";
    json_str += "};\n";

    json_str += "var trace_chh_body ={";
    json_str += "y:[" + list2string(chh_body_mr, chh_2) + "],";
    json_str += "name: 'chh_body',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#FF4136'}";
    json_str += "};\n";

    json_str += "var trace_chh_down ={";
    json_str += "y:[" + list2string(chh_down_mr, chh_1) + "],";
    json_str += "name: 'chh_down',";
    json_str += "type:'box',";
    json_str += "marker:{color:'#FF851B'}";
    json_str += "};\n";

    json_str += "var data = [trace_chh_up, trace_chh_body, trace_chh_down];\n";
    json_str += "var layout={title: 'CHH', xaxis: {tickangle: 45}, yaxis: { title:'Methylation level in CHH(%)', zeroline: false}, boxmode: 'group' };\n";

    json_str += "Plotly.newPlot('plot_enrich_chh', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    //enrich footer
    ofs << "</div>\n";
    ofs << "</div>\n";

    delete[] cg_down_mr;
    delete[] cg_body_mr;
    delete[] cg_up_mr;
    delete[] chg_down_mr;
    delete[] chg_body_mr;
    delete[] chg_up_mr;
    delete[] chh_down_mr;
    delete[] chh_body_mr;
    delete[] chh_up_mr;
    delete[] context;
    delete[] region;

    printFooter(ofs);
}

string list2string_long(long* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string list2string_int(int* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string list2string(double* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

// html function
void printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background: #21548e; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background: #21548e;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

void printJS(ofstream& ofs){
    //ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << endl;
    ofs << "<script src='https://www.dna-asmdb.com/download/plotly-latest.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string getCurrentSystemTime()
{
  //auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  //struct tm* ptm = localtime(&tt);
  char date[60] = {0};

  struct tm *current_date;
  time_t seconds;
  time(&seconds);
  current_date = localtime(&seconds);
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)1900+current_date->tm_year,(int)current_date->tm_mon+1,(int)current_date->tm_mday,
    (int)current_date->tm_hour,(int)current_date->tm_min,(int)current_date->tm_sec);
  return std::string(date);
}

void printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>batmeth2 report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "batmeth2 " << BT2_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
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

void plot_chrom(ofstream& ofs, int x_len, double* chrom_cg_p, double* chrom_cg_n,
    double* chrom_chg_p, double* chrom_chg_n, double* chrom_chh_p, double* chrom_chh_n, char* chrom){
    // script of plotly
    string chr = chrom;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "";
    const int xlen = x_len;
    int*  x_axis = new int[xlen];
    for(int i=0;i<x_len; i++){
        x_axis[i]=i;
    }

    json_str += "var trace_cg_p ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_cg_p, x_len) + "],";
    json_str += "name: 'strand: +',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#1C86EE', size: 2}";
    json_str += "};\n";

    json_str += "var trace_cg_n ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_cg_n, x_len) + "],";
    json_str += "name: 'strand: -',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#FF4136', size: 2}";
    json_str += "};\n";

    json_str += "var data = [trace_cg_p, trace_cg_n];\n";
    json_str += "var layout={title: 'CG', yaxis: { title:'DNA Methylation level(%)'}, xaxis: {title: ' " + chr + " '} };\n";

    json_str += "Plotly.newPlot('plot_chrom_cg_" + chr + "', data, layout);\n";

    //------non cg
    json_str += "var trace_chg_p ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_chg_p, x_len) + "],";
    json_str += "name: 'strand: +',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#1C86EE', size: 2}";
    json_str += "};\n";

    json_str += "var trace_chg_n ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_chg_n, x_len) + "],";
    json_str += "name: 'strand: -',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#FF4136', size: 2}";
    json_str += "};\n";

    json_str += "var data = [trace_chg_p, trace_chg_n];\n";
    json_str += "var layout={title: 'CHG', yaxis: { title:'DNA Methylation level(%)'}, xaxis: {title: ' " + chr + " '} };\n";

    json_str += "Plotly.newPlot('plot_chrom_chg_" + chr + "', data, layout);\n";

    //--- chh
    
    json_str += "var trace_chh_p ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_chh_p, x_len) + "],";
    json_str += "name: 'strand: +',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#1C86EE', size: 2}";
    json_str += "};\n";

    json_str += "var trace_chh_n ={";
    json_str += "x: [" + list2string_int(x_axis, x_len) + "],";
    json_str += "y: [" + list2string(chrom_chh_n, x_len) + "],";
    json_str += "name: 'strand: -',";
    json_str += "mode: 'markers',";
    json_str += "type:'scatter',";
    json_str += "marker: {color: '#FF4136', size: 2}";
    json_str += "};\n";

    json_str += "var data = [trace_chh_p, trace_chh_n];\n";
    json_str += "var layout={title: 'CHH', yaxis: { title:'DNA Methylation level(%)'}, xaxis: {title: ' " + chr + " '} };\n";

    json_str += "Plotly.newPlot('plot_chrom_chh_" + chr + "', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    
    delete []x_axis;
}

string my_to_string(long n)
{
    long m = n;
    int max=100;
    char s[max];
    char ss[max];
    int i=0,j=0;
    if (n < 0)// 处理负数
    {
        m = 0 - m;
        j = 1;
        ss[0] = '-';
    }    
    while (m>0)
    {
        s[i++] = m % 10 + '0';
        m /= 10;
    }
    s[i] = '\0';
    i = i - 1;
    while (i >= 0)
    {
        ss[j++] = s[i--];
    }    
    ss[j] = '\0';    
    return ss;
}

void outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + my_to_string(v) + "</td></tr>\n";
}

void outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}
