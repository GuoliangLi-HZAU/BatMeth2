BatMeth2: An Integrated Package for Bisulfite DNA Methylation Data Analysis with Indel-sensitive Mapping.  
--------------------------------------------------

This is a README file for the usage of Batmeth2.
--------------------------------------------------

Please go to https://github.com/GuoliangLi-HZAU/BatMeth2

1) batmeth2-master.zip (the zip of the program)

REQUIREMENTS
-------
1) gcc (v4.8) , gsl library, zlib < 1.2.5

2) R (ggplot2, pheatmap, xtable)

3) samtools (suggest: v1.3.1)

4) fastp, raw reads as input need

INSTALL
-------

a) Download 1) 

b) unzip 1) 

c) Change directory into the top directory of b) "BatMeth2/" 

d) Type 

- ./configure
- make
- make copy
- Rscript scripts/install.rpackages.r (Optional, used in visulization)


If your zlib version bigger than 1.2.5, and you do not need process file of gzip format, you can install the tools:
- ./configure
- make nogzip
- make copy-nogzip
- Rscript scripts/install.rpackages.r (Optional, used in visulization)


e) The binary of BatMeth2 will be created in bin/


BUILDING INDEX
-------

a) Have a fasta-formatted reference file ready 

b) Type "`BatMeth2 build_index GENOME.fa`"  for WGBS or `BatMeth2 build_index rrbs GENOME.fa` for RRBS to make the neccessary pairing data-structure based on FM-index.

c) Run "`BatMeth2`" to see information on usage.

USAGE of BatMeth2
------

**Example Data**

You can download the test data on `https://drive.google.com/open?id=1SEpvJbkjwndYcpkd39T11lrBytEq_MaC` .

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
BatMeth2 pipel --fastp ~/location/to/fastp -1 Raw_reads_1.fq.gz -2 Raw_read_2.fq.gz -g ./batmeth2index/genome.fa -o meth -p 6 --gff ./gene.gff -f 1
```

Or clean reads:

```bash
BatMeth2 pipel -1 Clean_reads_1.fq.gz -2 Clean_read_2.fq.gz -g ./batmeth2index/genome.fa -o meth -p 6 --gff ./gene.gff -f 1 
```


### BatMeth2 pipeline parameters

**BatMeth2 [mode][paramaters]**<br>

mode:  build_index, pipel, align, calmeth, annoation, methyPlot, batDMR, visul2sample, mkreport<br>

**[build_index]** <br>
    Usage:  (must run this step first)

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
    **[align paramaters:]** <br>
      see the details in 'BatMeth2 align' <br>
    **[calmeth paramaters:]** <br>
      see the details in 'BatMeth2 calmeth' <br>
    **[annotion paramaters:]** <br>
      see the details in 'BatMeth2 annoation' <br>
    **[methyPlot paramaters:]** <br>
      see the details in 'BatMeth2 methyPlot' <br>
    **[batDMR paramaters:]** <br>
      see the details in 'BatMeth2 batDMR' <br>
    **[visul2sample paramaters:]** <br>
      see the details in 'BatMeth2 visul2sample' <br>
    -h|--help   usage <br>


## OUTPUT FILE

Output file format and details see "https://github.com/GuoliangLi-HZAU/BatMeth2/blob/master/output_details.pdf".<br>
Output report details see "http://htmlpreview.github.io/?https://github.com/GuoliangLi-HZAU/BatMeth2/blob/master/BatMeth2-Report/batmeth2.html" or "http://htmlpreview.github.com/?https://github.com/GuoliangLi-HZAU/BatMeth2/blob/master/BatMeth2-Report/batmeth2.html" .<br>



## Functions in BatMeth2

BatMeth2 has the following main features: 

1. Batmeth2 has efficient and accurate alignment performance. 
2. Batmeth2 can calculate DNA methylation level of base site
3. BatMeth2 also can caculate and annotation DNA methylation level on chromosome region or gene/TE etc. functional region. 
4. By integrating BS-Seq data visualization (DNA methylation distribution on chromosome and gene etc) and BatMeth2 can show the results of the DNA methylation data more clearly. 
5. BatMeth2 can perform effective DNA methylation differential regions analysis based on the number of input samples and user requirements. And BatMeth2 provide differential methylation annotation ability.



#### 1. BatMeth2 alignment 

**Program: BatMeth2  align**

Version: v1.0

[ Single-end-reads ]

**Usage:**     batmeth2 -g INDEX -i INPUT -o OUTPUT

Example:   batmeth2 -g /data/index/hg19/hg19.fa -i Read.fq -o outPrefix -p 6 -n 2

[ Paired-end-reads ]

**Usage:**     batmeth2 -g INDEX -i INPUT_left -i INPUT_right -o OUTPUT

Example:   batmeth2 -g /data/index/hg19/hg19.fa -i Read_R1_left.fq -i Read_R2_right.fq -o outPrefix -p 6 -n 2

Parameters : 

​      --inputfile | -i <filename>   Name of input file

​      --genome | -g <filename>      Name of the genome mapped against

​      --outputfile | -o <filename>  Name of output file prefix

​      --buffersize | -b <integer>   Size of disk buffers

​      --indelsize                   indel size

​      --seed                        seed length for alignment. suggest: n/2. n = length of the read.

​      --maxhits | -m                maximum number of pairings to find

​      --maxmismatches | -n          maximum mismatches allowed due to seq. errors

​      --non_directional             Alignments to all four bisulfite strands will be reported. Default: OFF. 

​      --single | -u                 output singly mapped hits to file

​       --editdist | -e               maximum edit distance/percentage  to allow in mate rescue

​      --insertsize | -s             inital insert size, default 800

​      --matchlen | -E               minimum match len to allow in mate rescue

​      --maxtags | -t <integer>      maximum pairs to scan

​      --std | -d <integer>          standard deviatiion of reads distribution

​      --swon | -w(=rescdisc)        perform Smith-Waterman mate rescue. if rescdisc is set, try to rescue discordant pairs before marking them as discordant

​      --swlow | -Z                  Assign low quality to rescued mate

​      --flanksize | -f <integer>    size of flanking region for Smith-Waterman

​      --forcelength | -F <integer>  map only the first <integer> bases

​      --bwapol | -B                 descend only one mismatch mote than opt hit

​      --unique | -U                 Mark unique hits

​      --lowestmis | -P              Pair Hits having lowest count of mismatches in Head and Tail

​      --config | -c <file>          Load settings from configuration <file>  

​      --verbose | -v <number>       1-turn off general info 2-turn off progress bar 

​      --ntorand | N                 Replace N's in the read with random characters

​      --ncount | C                  Maximum N's allowed in the read

​      --logfile | -L                Log exit status to this log 

​      --heuristic                   Skip 2 mismatch scan of end anchors

​      --no_phred_filter             do not apply phred filter

​      --easyphred                   apply lenient phred filter

​      --general                     try to get most accurate mapping

​      --swlimit | <integer>         try at most <integer> sw extensions

​      --threads | -p <integer>      Launch <integer> threads

​      --NoInDels | -I               not to find the indels result

​      --help | -h                   Print help

Note: To use BatMeth2, you need to first index the genome with `build_all genome'.





#### 2. BatMeth2 calmeth

Command Format :  calmeth [options] -g GENOME  -i/-b <Samfile/Bamfile> -m <methratio outfile prefix> -p 6

**Usage:**

​       -g|--genome           Genome

​	-i|--input            Sam format file

​	-b|--binput           Bam format file

​	-p|--threads          the number of threads.

​	-n|--Nmismatch        Number of mismatches

​	-m|--methratio        [MethFileNamePrefix]  Predix of methratio output file

​	-Q [int]              caculate the methratio while read QulityScore >= Q. default:10

​	-c|--coverage         >= <INT> coverage. default:5

​	-nC		         >= <INT> nCs per region. default:5

​	-R |--Regions         Bins for DMR caculate , default 1kb .

​	--binsfile            DNA methylation level distributions in chrosome, default output file: {methratioPrefix}.methBins.txt

​	-s|--step             Chrosome using an overlapping sliding window of 100000bp at a step of 50000bp. default step: 50000(bp)

​	-r|--remove_dup       REMOVE_DUP, default:false

​	-f|--sam [outfile]    f for sam format outfile contain methState. default: sam format.

​	--sam-seq-beforeBS    Converting BS read to the genome sequences.

​	-h|--help 



#### 3. BatMeth2 annotation

BatMeth2:  MethyGff v1.0

Command Format :   methyGff [options] -o <OUT_PREFIX> -G GENOME -gff <GFF file>/-gtf <GTF file>/-b <bed file> -m <from Split methratio outfile> [-B][-P]

**Usage:**

​	-o|--out         Output file prefix

​	-G|--genome      Genome

​	-m|--methratio   Methratio output file.

​	-c|--coverage    >= <INT> coverage. default:5

​	-C               <= <INT> coverage. default 600.

​	-nC              >= <INT> Cs per bins or genes. default:5

​	-gtf|-gff        Gtf/gff file

​	-b|--BED         Bed file, chrom start end (strand)

​	-d|--distance    DNA methylation level distributions in body and <INT>-bp flanking sequences. The distance of upstream and downstream. default:2000

​	-B|--body        For different analysis input format, gene/TEs body methylation level. [Different Methylation Gene(DMG/DMT...)]

​	-P|--promoter    For different analysis input format.[Different Methylation Promoter(DMP)]

​	--TSS            Caculate heatmap for TSS. [Outfile: outPrefix.TSS.cg.n.txt]

​	--TTS            Caculate heatmap for TTS. [Outfile: outPrefix.TTS.cg.n.txt] 

​	--GENE           Caculate heatmap for GENE and flank 2k. [Outfile: outPrefix.GENE.cg.n.txt] 

​	-s|--step        Gene body and their flanking sequences using an overlapping sliding window of 5% of the sequence length at a step of 2.5% of the sequence length. So default step: 0.025 (2.5%)

​	-S|--chromStep   Caculate the density of genes/TEs in chromsome using an overlapping sliding window of 100000bp at a step of 50000bp, must equal "-s" in Split.. default step: 50000(bp)

​	-h|--help



#### 4. BatMeth2 methyPlot

**Usage1:**

* methyPlot chromsome.bins.txt chrosome.methy.distri.pdf step<default:0.025> Infile1.from.batmeth2:methyGff out1.pdf starLabel endLabel Infile2 out2.pdf

eg: methyPlot chromsome.bins.txt chrosome.methy.distri.pdf 0.025 gene.meth.Methylevel.1.txt methlevel.pdf TSS TTS gene.meth.AverMethylevel.1.txt elements.pdf

**Usage2: **

* methyPlot chromsome.bins.txt chrosome.methy.distri.pdf 0.025 gene.meth.Methylevel.1.txt methlevel.pdf TSS TTS gene.meth.AverMethylevel.1.txt elements.pdf test.annoDensity.1.txt test.density.pdf sampleElmentName test.mCdensity.txt test.mCdensity.pdf test.mCcatero.txt test.mCcatero.pdf 0.8 0.1 0.1



**Visulization case 2: two more samples**

*Contains:*

* The density of gene, transposon elements (TE) and the level of DNA methylation in the different samples of the whole genome. <br>
    $ `Rscript density_plot_with_methyl.r inputFile1 input2 genedensityFile TEdensity output.pdf label1 label2` 

    *example: Rscript ~/software/batmeth2/src/density_plot_with_methyl.r WT.methChrom.Rd.txt Mutant.methChrom.Rd.txt WT.noRd.Gff.gffDensity.1.txt WT.noRd.TE.gffDensity.1.txt density.Out.pdf WT mutant *

* DNA methylation level distribution across genes/TEs in different samples. 

    $ `Rscript methylevel.elements.r step(default:0.025) Input.Sample1.from.Batmeth2:methyGff Input.Sample2 outfilePrefix xLab1 xLab2 Sample1Prefix Sample2Prefix` 

    *example: Rscript methylevel.elements.compare.r 0.025 sample1.gene.meth.Methylevel.1.txt sample2.gene.meth.Methylevel.1.txt methlevel TSS TTS mutant WT * 



#### 5. BatMeth2 DMC or DMR/DMG

##### BatMeth2: DMR 

Command Format :  DMR [options] -g genome.fa -o_dm <DM_result>  -1 [Sample1-methy ..] -2 [sample2-methy ..] 

**Usage:**

​	-o_dm        output file

​	-o_dmr       when use auto detect by dmc

​	-g|--genome  Genome

​	-1           sample1 methy files, sperate by space.

​	-2           sample2 methy files, sperate by space.

​	-FDR         adjust pvalue cutoff default : 0.05

​	-methdiff    the cutoff of methylation differention. default: 0.25 [CpG]

​	-element     caculate gene or TE etc function elements.

​	-L          predefinded regions 

​	-h|--help



1. Pre-definded regions (Gene/TE/UTR/CDS...,but must run 'combined.element sample1 sample2 sample1out sample2out' before batDMR.) 

$ `BatMeth2 batDMR -g genome -L -o_dm dm.output.txt -1 [sample1.methC.txt replicates ..] -2 [sample2.methC.txt replicates ..]` 

2. Auto define DMR region according the dmc 

$ `BatMeth2 batDMR -g genome -o_dm dm.output.txt -o_dmr dmr.output.txt -1 [sample1.methC.txt replicates ..] -2 [sample2.methC.txt replicates ..]` 



##### The output format

```
    1. dmc <br> Chrom position starnd context pvalue adjust_pvalue combine_pvalue corrected_pvalue cover_sample1 meth_sample1 cover_sample2 cover_sample2 meth.diff 
    
    2. dmr <br> Chrom start end dmr score meth.diff aver_corrected_pvalue 
```

##### Filter the result

```bash
awk '$6<0.05 && sqrt($11*$11)>0.6 ' H9vsIMR90.gene.dmr.txt > H9vsIMR90.gene.dmr.filter.txt
```

##### DM annotation

BatMeth2:  DMCplot

Command Format :   DMCannotationPlot [options] -o <Out_File> -G GENOME -g <GFF files..> -d <dmc file> -c<mC context>

Usage:

​	-o|--out       Output file name.

​	-G|--genome    Genome file

​	-d|--dmcFile   dmc file. Format: Chrome Location strand

​                 Format: chr pos strand 

​	-g|--gff       Gff files, 1 or more

​	-c|--context   mC context, CG[default]/CHG/CHH/C. 

​	-h|--help 



```bash
DMCannotationPlot [options] -o <Out_File> -G GENOME -g <GFF files.. eg: TE.gff gene.gff CDS.gff intron.gff lncRNA.gff ...> -d <dmc file> -c <mC context default: CG>
```

   **Attention:**<br>
    *1.DMC file format : chr pos strand <br>
    2.GFF files are separated by spaces*<br><br>

------



Make sure all index files reside in the same directory.

Built with `BatMeth2 build_index Genome.fa` 

=-=-=-=-=-=-=-=-=-=

GNU automake v1.11.1, GNU autoconf v2.63, gcc v4.4.7.

Tested on Red Hat 4.4.7-11 Linux

Thank you for your patience.
