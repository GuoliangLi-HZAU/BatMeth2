#ifndef _SAM_FUNCS_H
#define _SAM_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <bam.h>

#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4
#endif

#ifndef BAM_CIGAR_MASK
#define BAM_CIGAR_MASK 0xf
#endif

#include "hash_funcs.h"

typedef struct MapRecord_Struct MapRecord;
struct MapRecord_Struct
{
	int readId;		// Not used
	char *qname;
	char strand;
	char *chrome;
	int flag;
	int offset;
	unsigned int mapq;
	char *cigar;
	char *seq;
	char *seqBuf;
	char *qual;
	char *qualBuf;
	char *comBuf;
	int len;
	int r12;
};

void snpAnalysis(char* bamFileName, char* snpFileName, char* methCgFileName, char* methChgFileName, char* methChhFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int vQualMin, int nLayerMin, int nLayerMax, float vSnpRate, float vSnpPerBase, unsigned int mapqThr);
void methProcess(char* bamFileName, char* methCgFileName, char* methChgFileName, char* methChhFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int vQualMin, int nLayerMin, unsigned int mapqThr);
void snpProcess(char* bamFileName, char* snpFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int chrCnt, int vQualMin, int nLayerMax, float vSnpRate, float vSnpPerBase, unsigned int mapqThr);

int parseBuffer(bam_header_t *header, bam1_t *b, MapRecord* record, unsigned int mapqThr);
char nxtChar(char* cigar, int* len);
void dispRecord(MapRecord* record);
int seqReverse(char* seq, int seqLen);
void printMeth(FILE* methFptr, int len, char* curChr, unsigned int* w_Mm, unsigned int* w_Mc, unsigned int* w_Mq, unsigned int* c_Mm, unsigned int* c_Mc, unsigned int* c_Mq, char* tag, int nLayerMin);
void printSnp(FILE* posFptr, char** chrSeqArray, int idx, int len, float vSnpRate, char* curChr, unsigned short *w_A, unsigned short *w_T, unsigned short *w_C, unsigned short *w_G, unsigned short *c_A, unsigned short *c_T, unsigned short *c_C, unsigned short *c_G, unsigned int *w_Aq, unsigned int *w_Tq, unsigned int *w_Cq, unsigned int *w_Gq, unsigned int *c_Aq, unsigned int *c_Tq, unsigned int *c_Cq, unsigned int *c_Gq, unsigned short *w_An, unsigned short *w_Tn, unsigned short *w_Cn, unsigned short *w_Gn, unsigned short *c_An, unsigned short *c_Tn, unsigned short *c_Cn, unsigned short *c_Gn, unsigned int *w_Q, unsigned int *c_Q);

#endif
