"BatMeth2: ...".  
--------------------------------------------------

This is a README file for the usage of Batmeth2.
--------------------------------------------------

Please go to https://github.com/ZhouQiangwei/batmeth2

1) batmeth2-master.zip (the zip of the program)
   
   
INSTALL
-------
   
a) Download 1) 

b) unzip 1) 

c) Change directory into the top directory of b) "batmeth2/" 

d) Type " ./configure " then type " make " and finally " make copy " 

e) The binary of Batmeth2 will be created in bin/
   

BUILDING INDEX
-------
   
a) Have a fasta-formatted file ready 

b) Locate the script "build_all" in "batmeth2/bin" 

c) Type "`build_all GENOME.fa`" to make the neccessary pairing data-structure based on FM-index.
   

USAGE
------
   
**[1] Mapping** 

**Single-end-reads** 

Short reads align (length < 150) 

$ `batmeth2-align -g INDEX -i INPUT -o OUTPUT_Prefix -p 6 -n 2` 

*example: 
    1. Left: batmeth2-align -g ./hg19/hg19.fa -i left.fq -o outPrefix1 -p 6 <br>
    2. Right: batmeth2-align --non_directional -O -g ./hg19/hg19.fa -i right.fq -o outPrefix2 -p 6 
*

Long reads align (length > 150)

$ ` batmeth2-l -g INDEX -i INPUT -o OUTPUT_Prefix -p 6 -n 2 ` 


**Paired-end-reads** 

$ `batmeth2-align -g INDEX -i INPUT_left -i INPUT_right -o OUTPUT_prefix --threads 16 -n 2` 

*example: batmeth2-align -g /data/index/hg19/hg19.fa -i CML_R1_left.fq -i CML_R2_right.fq -o out.prefix --threads 6 
*
   

**[2] Get Methylation**

$ `split [options] -o [FINAL_RESULT] -g GENOME -n [number of mismatches] -i [all temp(s) from batmeth] -m [methratio outfile]`

example: split -f -o Final_Result.sam -g ../../../Genome/batmeth2/all.con -n 2 -i batmeth2outPrefix.* -m Final.methratio.txt
   
**[3] Detect BSseq SNP using approximate Bayesian modeling** 

a) Convert Sam files into BAM files

   $ `samtools view -bS Final_Result.sam > Final_Result.bam` 
   

b) Sort the bam file 

   $ `samtools sort Final_Result.bam -f Final_Result.sort.bam` 
   
   $ `rm Final_Result.bam` 
   
c) BS-Seq variation detection (use BS-Snper software)

   $ `BS-Snper --fa GENOME --input Sort.bam --output snp_result_file` 
   
   *example: BS-Snper --fa /data/index/hg19/hg19.fa --input Final_Result.sort.bam --output snp_result_file* 
   
**[4] Gene OR TEs DNA methylation level and density** 

   $ `methyGff [options] -o [OUT_PREFIX] -G GENOME -g <GFF file>/-b <bed file> -m [from Split methratio outfile] [-B] [-P]` 
   
**[5] DNA methylation data visulization**  *(required R package : ggplot2{install.packages("ggplot2")})*  

   1. Only one sample. 

    *Contains:* 
    * DNA methylation density level of chromsome.  
    
    * DNA methylation distribution across gene/TE. 
    
    * The DNA methylation level in 2Kb up of gene, gene body, 2Kb down of gene body. 
   
        $ `methyPlot chromsome.bins.txt chrosome.methy.distri.pdf step[default:0.025] Infile1.from.batmeth2:methyGff out1.pdf Infile2 out2.pdf ` 
        
        *example: methyPlot chromsome.bins.txt chrosome.methy.distri.pdf 0.025 gene.meth.Methylevel.1.txt methlevel.pdf gene.meth.AverMethylevel.1.txt elements.pdf* 
   
    * The density of gene, transposon elements (TE) and the level of DNA methylation in the whole genome. 
    
        $ `Rscript density_plot_with_methyl.r inputFile1 genedensityFile TEdensity output.pdf label` 
        
        *example: Rscript ~/software/batmeth2/src/density_plot_with_methyl_oneSample.r methChrom.Rd.txt H9.noRd.Gff.gffDensity.1.txt H9.noRd.TE.gffDensity.1.txt density.test.pdf H9no* 

    * DNA methylation heatmap in all genes. 
    
        $ `GeneMethHeatmap Wt[the prefix of methyGff output] mutant CG-ceil CHG-ceil CHH-ceil` 
        
        *example: GeneMethHeatmap H9.Rd.Gff IMR90.Rd.Gff 1.0 0.6 0.2* 
        

   2. Two more samples. 

    *Contains:*
    * The density of gene, transposon elements (TE) and the level of DNA methylation in the different samples of the whole genome. <br>
        $ `Rscript density_plot_with_methyl.r inputFile1 input2 genedensityFile TEdensity output.pdf label1 label2` 
        
        *example: Rscript ~/software/batmeth2/src/density_plot_with_methyl.r WT.methChrom.Rd.txt Mutant.methChrom.Rd.txt WT.noRd.Gff.gffDensity.1.txt WT.noRd.TE.gffDensity.1.txt density.Out.pdf WT mutant *
    * DNA methylation level distribution across genes/TEs in different samples. 
    
        $ `Rscript methylevel.elements.r step(default:0.025) Input.Sample1.from.Batmeth2:methyGff Input.Sample2 outfilePrefix xLab1 xLab2 Sample1Prefix Sample2Prefix` 
        
        *example: Rscript methylevel.elements.compare.r 0.025 sample1.gene.meth.Methylevel.1.txt sample2.gene.meth.Methylevel.1.txt methlevel TSS TTS mutant WT * 
        

**[6] Get DMC or DMR/DMG**  

   1. Pre-definded regions (Gene/TE/UTR/CDS...,but must run 'combined.element sample1 sample2 sample1out sample2out' before batDMR.) 

        $ `batDMR -g genome -L -o_dm dm.output.txt -1 [sample1.methC.txt replicates ..] -2 [sample2.methC.txt replicates ..]` 
        
        *example: batDMR -g genome.fa -L -o_dm H9vsIMR90.gene.dmr.txt -1 H9.gene.meth.txt -2 IMR90.gene.meth.txt *  
        
   2. Auto define DMR region according the dmc 

        $ `batDMR -g genome -o_dm dm.output.txt -o_dmr dmr.output.txt -1 [sample1.methC.txt replicates ..] -2 [sample2.methC.txt replicates ..]` 
        
        *example: batDMR -g genome.fa -o_dm H9vsIMR90.dmc.txt -o_dmr H9vsIMR90.dmr.txt -1 H9.chr10.methC.txt -2 IMR90.chr10.methC.txt *  
        
    *The output format:* 
    
        *1. dmc <br> Chrom position starnd context pvalue adjust_pvalue combine_pvalue corrected_pvalue cover_sample1 meth_sample1 cover_sample2 cover_sample2 meth.diff * 
        
        *2. dmr <br> Chrom start end dmr score meth.diff aver_corrected_pvalue * 
        
    *filter the result:* 
    
        $ `awk '$6<0.05 && sqrt($13*$13)>0.6 ' H9vsIMR90.gene.dmr.txt > H9vsIMR90.gene.dmr.filter.txt` 
        

**[7] DM annotation** 

   $ `DMCannotationPlot [options] -o <Out_File> -G GENOME -g <GFF files.. eg: TE.gff gene.gff CDS.gff intron.gff lncRNA.gff ...> -d <dmc file> -c <mC context default: CG>` 
   
   **Attention:**<br>
    *1.DMC file format : chr pos strand <br>
    2.GFF files are separated by spaces*<br><br>

INDEX was the mentioned "hg19.fa". Make sure all index files reside in the same directory.

Built with `build_all Genome.fa` 

=-=-=-=-=-=-=-=-=-=

GNU automake v1.11.1, GNU autoconf v2.63, gcc v4.4.7.

Tested on Red Hat 4.4.7-11 Linux

Thank you for your patience.
