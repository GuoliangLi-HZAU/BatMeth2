/* prints bases and quals */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bam.h>

#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4
#endif

#ifndef BAM_CIGAR_MASK
#define BAM_CIGAR_MASK 0xf
#endif

int main(int argc,char** argv)
{
	int i;
	uint8_t *s, *t;
    bam_header_t *header;
	
	bamFile in = bam_open(argv[1], "r");
    if(in == NULL) {
		fprintf(stderr, "Cannot open bam file!\n");
		return -1;
	}
	
    bam1_t* b=bam_init1();
    if(b == NULL) {
		fprintf(stderr, "Cannot init bam structure!\n");
		return -1;
	}
	
    header = bam_header_read(in);
	while(bam_read1(in, b) >= 0) {
		const bam1_core_t *c = &b->core;
		s = bam1_seq(b);
		t = bam1_qual(b);


		// 1 - QNAME
		fprintf(stderr, "%s\t", bam1_qname(b));
		// 2 - FLAG
		fprintf(stderr, "%d\t", c->flag);
		// 3 - RNAME
		fprintf(stderr, "%s\t", header->target_name[c->tid]);
		// 4 - POS
		fprintf(stderr, "%d\t", c->pos);
		// 5 - MAPQ
		fprintf(stderr, "%d\t", (int)(c->qual));
		// 6 - CIGAR
		for (i = 0; i < c->n_cigar; ++i)
			fprintf(stderr, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
		fprintf(stderr, "\t");
		
		// 10 - SEQ
		for(i = 0; i < c->l_qseq; ++i) 
			fprintf(stderr, "%c", bam_nt16_rev_table[bam1_seqi(s, i)]);
		fprintf(stderr, "\t");
		
		// 11 - QUAL
        if(t[0] == 0xff) {
			fprintf(stderr, "*");
		}
		else {
			fprintf(stderr, "------------------------\n");
			for(i = 0; i < c->l_qseq; ++i) 
				fprintf(stderr, "%c\t", t[i] + 33);	
			fprintf(stderr, "\n");	
			for(i = 0; i < c->l_qseq; ++i) 
				fprintf(stderr, "%u\t", (unsigned short)t[i]);
			fprintf(stderr, "------------------------\n");
		}
		fprintf(stderr, "\n");
		
		return 0;
	}
	
	bam_header_destroy(header);
	bam_close(in);
	bam_destroy1(b);
	
	return 0;
 }