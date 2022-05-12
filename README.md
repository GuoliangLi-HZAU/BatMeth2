BatMeth2: An Integrated Package for Bisulfite DNA Methylation Data Analysis with Indel-sensitive Mapping.  
--------------------------------------------------

This is a README file for the usage of Batmeth2.
--------------------------------------------------

BatMeth2 tutotial: https://batmeth2-docs.readthedocs.io

Starting from this version, because the alignment problem needs to be modified for a long time, we replaced the align part written by ourselves with BWA MEM, and added some Python visualization functions. At present, some functions are being added.

REQUIREMENTS
-------
1) gcc (v4.8) , gsl library, zlib

2) samtools >= v1.3.1

3) fastp, raw reads as input need

INSTALL
-------

a) Download 1) 

b) unzip 1) 

c) Change directory into the top directory of b) "BatMeth2/" 

d) Type 

- ./configure
- make
- make install

e) The binary of BatMeth2 will be created in bin/


Building index
-------

a) Have a fasta-formatted reference file ready 

b) Type "`BatMeth2 build_index GENOME.fa`"  for WGBS or `BatMeth2 build_index rrbs GENOME.fa` for RRBS to make the neccessary pairing data-structure based on FM-index.

c) Run "`BatMeth2`" to see information on usage.

Usage of BatMeth2
------

**Example Data**

You can download the test data on `https://drive.google.com/open?id=1SEpvJbkjwndYcpkd39T11lrBytEq_MaC` 
Or `https://pan.baidu.com/s/1mliGjbn_33wlQLieqy5YOQ` with extraction code: `kr32`.

**Example data contain files:**

* input fastq.gz (paired end)
* genome file
* usage code and details
* gene annotation file

### Citation:
[Zhou Q, Lim J-Q, Sung W-K, Li G: An integrated package for bisulfite DNA methylation data analysis with Indel-sensitive mapping. BMC Bioinformatics 2019, 20:47.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2593-4)


### BatMeth2 pipeline

An easy-to-use, auto-run package for DNA methylation analyses

​    In order to complete the DNA methylation data analysis more conveniently, we packaged all the functions to complete an easy-to-use, auto-run package for DNA methylation analysis. During the execution of BatMeth2 Tool, an html report is generated about statistics of the sample.

The usage is here:

Raw reads:

```bash
BatMeth2 pipel --fastp ~/location/to/fastp -1 Raw_reads_1.fq.gz -2 Raw_read_2.fq.gz -g ./batmeth2index/genome.fa -o meth -p 8 --gff ./gene.gff
```

Or clean reads:

```bash
BatMeth2 pipel -1 Clean_reads_1.fq.gz -2 Clean_read_2.fq.gz -g ./batmeth2index/genome.fa -o meth -p 8 --gff ./gene.gff
```


### BatMeth2 pipeline parameters

**BatMeth2 [mode][paramaters]**

mode:  build_index, pipel, align, calmeth, annoation, batDMR<br>

**[build_index]** <br>
    Usage:  (must run this step first) <br>

1. BatMeth2 build_index genomefile.<br>

2. BatMeth2 build_index rrbs genomefile. <br>

**[pipel (Contains: align, calmeth, annoation, methyPlot, mkreport)]** <br>
    **[fastp location]** <br>
      &emsp;&emsp; --fastp    fastp program location. <br>
      &emsp;&emsp; If --fastp is not defined, the input file should be clean data. <br>
    **[select aligner]** <br>
      --aligner &emsp;&emsp; BatMeth2(default), bwa-meth, bsmap, bismark2, no (exit output_prefix.sam file, no need align again) <br>
    **[other aligners paramaters]** <br>
      --go &emsp;&emsp;&emsp; Name of the genome, contaion index build by aligner. (bwa-meth/bismark2) <br>
    **[main paramaters]** <br>
       --config &emsp; [config file]. &emsp;&emsp;  When we run pipel function in batches datasets, <br>
       &emsp;&emsp;&emsp;&emsp; please fill in the specified configuration file. <br>
       &emsp;&emsp;&emsp;&emsp; And there is a sample file (multirun.onf) in the BatMeth2 directory. <br>
       --mp&emsp;[4] &emsp; When batch processing data, we set the number of samples to run at a time (-mp, default is 4), and each sample needs six threads (- P parameter) by default. <br>
      -o &emsp; Name of output file prefix<br>
      -O &emsp; Output of result file to specified folder, default output to current folder (./) <br>
    **[alignment paramaters]** <br>
      -i &emsp;&ensp; Name of input file, if paired-end. please use -1, -2, input files can be separated by commas <br>
      -1 &emsp;&ensp; Name of input file left end, if single-end. please use -i <br>
      -2 &emsp;&ensp; Name of input file left end <br>
      -g &emsp; Name of the genome mapped against <br>
      -n &emsp; maximum mismatches allowed due to seq. errors <br>
      -p <interger> &emsp; Launch <integer> threads <br>
    **[calmeth paramaters]** <br>
      --Qual &emsp;&ensp;&ensp; calculate the methratio while read QulityScore >= Q. default:10 <br>
      --redup &emsp;&ensp; REMOVE_DUP, 0 or 1, default 0 <br>
      --region &emsp; Bins for DMR calculate , default 1000bp . <br>
      -f  &emsp;&ensp; for sam format outfile contain methState. [0 or 1], default: 0 (dont output this file). <br>
    **[calmeth and annoation paramaters]** <br>
      --coverage &emsp;&ensp; >= <INT> coverage. default:5 <br>
      --binCover &emsp;&ensp; >= <INT> nCs per region. default:3 <br>
      --chromstep &emsp; Chromosome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp) <br>
    **[annoation paramaters]** <br>
      --gtf/--gff/--bed &emsp; Gtf or gff file / bed file <br>
      --distance &emsp;&ensp; DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000 <br>
      --step &emsp;&emsp; Gene body and their flanking sequences using an overlapping sliding window of 5% of the sequence length at a step of 2.5% of the sequence length. So default step: 0.025 (2.5%) <br>
      -C   &emsp;&emsp;   <= <INT> coverage. default:1000 <br>
    **[mkreport paramaters]** <br>
      Make a batmeth2 html report, can see the detail in BatMeth2_Report/ directory. <br>
      **-o [outprefix]** <br>
    -h|--help   usage <br>


## Output files

Output file format and details see "https://github.com/GuoliangLi-HZAU/BatMeth2/blob/master/output_details.pdf".<br>

Output report details see "https://www.dna-asmdb.com/download/batmeth2.html" .<br>


## Functions in BatMeth2

BatMeth2 has the following main features: 

1. Batmeth2 has efficient and accurate alignment performance. 
2. Batmeth2 can calculate DNA methylation level of base site
3. BatMeth2 also can caculate and annotation DNA methylation level on chromosome region or gene/TE etc. functional region. 
4. By integrating BS-Seq data visualization (DNA methylation distribution on chromosome and gene etc) and BatMeth2 can show the results of the DNA methylation data more clearly. 
5. BatMeth2 can perform effective DNA methylation differential regions analysis based on the number of input samples and user requirements. And BatMeth2 provide differential methylation annotation ability.



---------------------

Make sure all index files reside in the same directory.

Built with `BatMeth2 build_index Genome.fa` 

=-=-=-=-=-=-=-=-=-=

GNU automake v1.11.1, GNU autoconf v2.63, gcc v4.4.7.

Tested on Red Hat 4.4.7-11 Linux

Thank you for your patience.
